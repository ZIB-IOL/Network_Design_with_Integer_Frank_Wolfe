function initialize_parameters(ta_data_name, alg, base_cost, rem_cd, total_scenarios, scale)
    # Load transportaiton Network
    tn = load_ta_network(ta_data_name) 

    if rem_cd == true
        @info "Removing circular demand"
        tn = remove_circular_demand(tn)
    end

    # The followign lines are to extract the information from the network and create the relevant data structures. 

    global MAX_FLOW = 1.1*tn.total_od_flow

    node_count = tn.number_of_nodes
    edge_count = tn.number_of_links

    if tn.first_thru_node > 1
        zone_count = tn.first_thru_node -1 
    elseif tn.first_thru_node == 1 && tn.number_of_zones == node_count
        zone_count = 1
    elseif tn.first_thru_node == 1 && tn.number_of_zones != node_count
        zone_count = tn.number_of_zones
    else
        throw(ArgumentError("Incorrect Argument for method. $(alg)"))
    end

    base_travel_demands = tn.travel_demand

    @assert all(base_travel_demands .>= 0)
    @assert any(base_travel_demands .> 0)
    @assert all(tn.capacity .> 0)
    @assert all(tn.power .>= 0)

    commodity_types = node_count

    od_pair_count = zone_count > 1 ? zone_count : node_count

    if total_scenarios == 1
        scenarios = Dict([s=>Dict([i => 1 .+ zeros(commodity_types) for i in range(1, step=1, stop=od_pair_count)]) for s in 1:total_scenarios])
    else
        scenarios = Dict([s=>Dict([i => 1 .+ scale * rand(commodity_types) for i in range(1, step=1, stop=od_pair_count)]) for s in 1:total_scenarios])
    end
    #CSV.write("scenarios_$(alg)_$(ta_data_name)_$(total_scenarios).csv", scenarios)

    #node_demands_dict = Dict([i=>zeros(commodity_types) for i in range(1,step=1,stop=od_pair_count)])
    node_demands_dict = Dict([s=>Dict([i => zeros(commodity_types) for i in range(1, step=1, stop=od_pair_count)]) for s in 1:total_scenarios])
    travel_demands = Dict([s =>  similar(tn.travel_demand) for s in 1:total_scenarios])

    for scenario_cnt in 1:total_scenarios
        for dest in range(1, step=1, stop=od_pair_count)
            for origin in range(1, step=1, stop=od_pair_count)
                #if origin == dest
                #    node_demands_dict[dest][origin] = -sum(travel_demands[:, dest])
                #else
                link_demand = base_travel_demands[origin, dest]
                node_demands_dict[scenario_cnt][dest][origin] = round(link_demand * scenarios[scenario_cnt][dest][origin], digits=3)
                travel_demands[scenario_cnt][origin, dest] = round(link_demand * scenarios[scenario_cnt][dest][origin], digits=3)
                #end
            end
            if travel_demands[scenario_cnt][dest, dest] > 0.0
                @info "There is circular demand for node $(dest)."
            end
        end
    end

    @info "Creating Graph"
    g = Graphs.SimpleDiGraph(node_count)

    edge_costs_dict = Dict{Tuple{Int64, Int64}, Vector{Float64}}()
    edge_capacities_dict = Dict{Tuple{Int64, Int64}, Vector{Float64}}()
    edge_cum_cap_dict = Dict{Tuple{Int64, Int64}, Float64}()

    g1 = MC_graph_with_weights(commodity_types, g, node_demands_dict, edge_costs_dict, edge_capacities_dict, edge_cum_cap_dict, zone_count)

    for edge in range(1, step=1, stop=edge_count)
        src = tn.init_node[edge]
        dst = tn.term_node[edge]

        weights = tn.link_length[edge] * ones(commodity_types)
        cap = tn.capacity[edge]
        add_mc_edge_wts!(g1, src, dst, commodity_types, weights, cap*ones(commodity_types), cap)
    end


    no_edges = Graphs.ne(g1.graph) 
    no_nodes = Graphs.nv(g1.graph)
    noz = od_pair_count
    edge_list = [(tn.init_node[i], tn.term_node[i]) for i in 1:Graphs.ne(g1.graph)]
    edge_dict = Dict([edge_list[i]=>i for i in eachindex(edge_list)])
    incoming_edges, outgoing_edges = neighbouring_edges(g1.graph)

    @info "Created graph with $(no_nodes) nodes, $(no_edges) edges and $(noz) zones"

    params = Params(no_edges, no_nodes, noz, 
                    edge_list, edge_dict, incoming_edges, 
                    outgoing_edges, node_demands_dict, MAX_FLOW)

    init_nodes = [i for i in tn.init_node]
    term_nodes = [j for j in tn.term_node]

    first_thru_node = tn.first_thru_node
    link_dic = sparse(init_nodes, term_nodes, collect(1:no_edges))

    gparams = GraphParams(
        init_nodes,
        term_nodes,
        travel_demands,
        first_thru_node,
        link_dic,
        od_pair_count,
        [],
        [],
        1.1 * tn.total_od_flow,
        total_scenarios
    )

    return tn, g1, params, gparams
    
end

"""
This function initimizes the optimizer and the constraints 
for the flow polytope.
"""
function initialize_model(alg, tn, g1, params, gparams, tfc, total_scenarios, solver=nothing)
    @info "Initializing the model according to the algorithm $(alg)."

    if (alg == "IFW" || alg == "LP" || alg == "SCIP" || alg == "IFW-P")
        #time_limit = 1800.0
        mip_rel_gap = 1e-4

        if solver == "SCIP"
            optimizer = SCIP.Optimizer(display_verblevel=0, limits_gap=mip_rel_gap)
        elseif solver == "Gurobi"
            optimizer = Gurobi.Optimizer()
            #MOI.set(optimizer, MOI.RawOptimizerAttribute("TimeLimit"), time_limit)
            MOI.set(optimizer, MOI.RawOptimizerAttribute("MIPGap"), mip_rel_gap)
        elseif solver == "HiGHS"
            optimizer = HiGHS.Optimizer()
            #MOI.set(optimizer, MOI.RawOptimizerAttribute("time_limit"), time_limit)
            MOI.set(optimizer, MOI.RawOptimizerAttribute("mip_rel_gap"), mip_rel_gap)
        else
            throw(ArgumentError("Invalid Solved Specified for IFW Method $(solver)"))
        end
        
        #optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 1800, "MIPGap" => 1e-3)
        #MOI.empty!(optimizer)
        

        if tn.number_of_zones == tn.number_of_nodes
            @info "The number of zones and nodes are the same."
            for scenario_cnt in 1:total_scenarios
                optimizer, con_list = create_flow_polytope(optimizer, g1, params, scenario_cnt)
            end
        else
            @info "The number of zones and nodes are different."
            for scenario_cnt in 1:total_scenarios
                optimizer, con_list = create_flow_polytope_zones(optimizer, g1, params, scenario_cnt)
            end
        end

    elseif (alg == "NLMO-IFW" || alg == "NLMO-P")

        optimizer = NDFW.Optimizer()


        optimizer.graph = g1.graph
        optimizer.true_graph = g1.graph
        optimizer.gparams = gparams
        optimizer.curr_binary_soln = ones(length(gparams.removed_edges))
        optimizer.type_feas_cuts = tfc
        optimizer.method = alg

        if tn.number_of_zones == tn.number_of_nodes
            @info "The number of zones and nodes are the same."
            for scenario_cnt in 1:total_scenarios
                optimizer.sopt, con_list = create_flow_polytope(optimizer.sopt, g1, params, scenario_cnt)
            end
        else
            @info "The number of zones and nodes are different."
            for scenario_cnt in 1:total_scenarios
                optimizer.sopt, con_list = create_flow_polytope_zones(optimizer.sopt, g1, params, scenario_cnt)
            end
        end

    elseif alg == "NLMO-FW" || alg == "BNDLMO-P" || alg == "BDM"
        optimizer = nothing
        con_list = Dict()
    else
        throw(ArgumentError("Incorrect Argument for method. $(alg)"))
    end

    return optimizer, con_list
end


function initialize_lmo(alg, optimizer, g, gparams, use_bigm, use_adaptive, use_reverse)
    if (alg == "IFW" || alg == "FW" || alg == "IFW-P")
        lmo = Boscia.MathOptBLMO(optimizer, false)
    elseif alg == "NLMO-FW"
        lmo = NetworkLMO(g, gparams)
    elseif alg == "NLMO-IFW"
        lmo = NDFW.MOINetworkLMO(optimizer)
    elseif alg == "NLMO-P"
        lmo = NDFW.NetworkLMO_P(optimizer)
    elseif alg == "BNDLMO-P"
        lmo = initialize_bounded_lmo(g, gparams, use_bigm, use_adaptive, use_reverse)
    elseif alg == "BDM"
        lmo = initialize_bdlmo(g, gparams, use_adaptive)
    else
        throw(ArgumentError("Incorrect Argument for method. $(alg)"))
    end

    return lmo
end

"""
Initialize the bounded LMO with the graph, the bounds and the integer variables
"""
function initialize_bounded_lmo(g, gparams, use_bigm, useadaptive, use_reverse)
    bounds = Boscia.IntegerBounds()

    od_pair_count = gparams.od_pair_count
    edge_count = length(gparams.init_nodes)
    new_edge_count = length(gparams.removed_edges)
    total_scenarios = gparams.total_scenarios

    total_scenario_vars = od_pair_count * edge_count + edge_count

    no_vars = total_scenario_vars * total_scenarios + new_edge_count
    # for i in 1:(total_scenario_vars*total_scenarios)
    #     push!(bounds, (i, 0.0), :greaterthan)
    #     push!(bounds, (i, Inf), :lessthan)
    # end

    int_vars = collect((total_scenario_vars*total_scenarios+1):no_vars)

    for i in int_vars
        push!(bounds, (i, 0.0), :greaterthan)
        push!(bounds, (i, 1.0), :lessthan)
    end

    lmo = NDFW.BNDesignLMO(no_vars, g, int_vars, bounds, gparams, use_bigm, useadaptive, use_reverse)

    return lmo
end


"""
Initialize the bounded LMO with the graph, the bounds and the integer variables
"""
function initialize_bdlmo(g, gparams, useadaptive)
    bounds = Boscia.IntegerBounds()

    od_pair_count = gparams.od_pair_count
    edge_count = length(gparams.init_nodes)
    new_edge_count = length(gparams.removed_edges)
    total_scenarios = gparams.total_scenarios
    nc = Graphs.nv(g)
    distmx = Inf * ones(nc, nc)

    total_scenario_vars = od_pair_count * edge_count + edge_count

    no_vars = total_scenario_vars * total_scenarios + new_edge_count
    # for i in 1:(total_scenario_vars*total_scenarios)
    #     push!(bounds, (i, 0.0), :greaterthan)
    #     push!(bounds, (i, Inf), :lessthan)
    # end

    int_vars = collect((total_scenario_vars*total_scenarios+1):no_vars)

    for i in int_vars
        push!(bounds, (i, 0.0), :greaterthan)
        push!(bounds, (i, 1.0), :lessthan)
    end

    lmo = NDFW.BDDesignLMO(no_vars, g, int_vars, bounds, gparams, useadaptive, distmx)

    return lmo
end