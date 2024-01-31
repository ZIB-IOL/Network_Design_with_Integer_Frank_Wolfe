"""
This function initimizes the optimizer and the constraints 
for the flow polytope.
"""
function initialize_stochastic_model(alg, tn, g1, params, gparams, tfc)
    @info "Initializing the model according to the algorithm $(alg)."

    if (alg == "IFW" || alg == "FW" || alg == "SCIP" || alg == "IFW-P")
        optimizer = SCIP.Optimizer(limits_time=800, display_verblevel=1, limits_gap=1e-4)

        if tn.number_of_zones == tn.number_of_nodes
            @info "The number of zones and nodes are the same."
            optimizer, con_list = create_flow_polytope(optimizer, g1, params)
        else
            @info "The number of zones and nodes are different."
            optimizer, con_list = create_flow_polytope_zones(optimizer, g1, params)
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
            optimizer.sopt, con_list = create_flow_polytope(optimizer.sopt, g1, params)
        else
            @info "The number of zones and nodes are different."
            optimizer.sopt, con_list = create_flow_polytope_zones(optimizer.sopt, g1, params)
        end

    elseif alg == "NLMO-FW"
        optimizer = nothing
        con_list = []
    else
        throw(ArgumentError("Incorrect Argument for method. $(alg)"))
    end

    return optimizer, con_list
end

"""
This function runs the Traffic assignment problem with the Frank Wolfe algorithm 
with and without the network design constraints. When network design constraints 
are present it uses Boscia.jl and when they are not present it uses blended_pairwise_conditional_gradient from the 
FrankWolfe.jl library. Note that the arcs which are to be expanded are determined at random. 
INPUT
ta_data_name : Name of data set from the TransportationNetworks library
enable_network_design : whether to use the network design constraints
base_cost : cost per link of network expansion
res: solution of the tailored frank wolfe algorithm. to be used for scaling. 
OUTPUT:
The above three inputs ans a result which constains the optimal objective function
"""
function stochastic_network_design(cost_of_expansion, sres, econfig, scenarios, scenario_count, edges2remove = 0)
    ta_data_name = econfig.ta_data_name
    alg = econfig.alg
    base_cost = cost_of_expansion
    res = sres
    time_limit = econfig.time_limit
    max_fw_iter = econfig.max_fw_iter
    remove_circular_demand = econfig.remove_circular_demand
    type_feas_cuts = econfig.type_feas_cuts
    fraction_removed = econfig.fraction_removed
    bs_verb = econfig.bs_verb
    fw_verb = econfig.fw_verb
    ls_verb = econfig.ls_verb
    dual_tightening = econfig.dual_tightening
    rel_dual_gap = econfig.rel_dual_gap

    tn, g1, params, gparams = initialize_parameters(ta_data_name, alg, base_cost, remove_circular_demand)

    #The following lines use the previously defined parameter values to crate the optimization problem. 


    optimizer, con_list = initialize_model(alg, tn, g1, params, gparams, type_feas_cuts)
    
    cost_of_expansion = []


    if (alg == "IFW" || alg == "NLMO-IFW")
        # chose 50% of arcs on which allow for network expansion
        seed = econfig.seed
        @info "Choosing random seed $(seed)"
        Random.seed!(seed)
        
        if edges2remove == 0
            #arcs_to_be_expanded = [rand() > 0.5 for i in 1:length(edge_list)]
            count_arcs_to_remove = Int(ceil(fraction_removed * length(params.edge_list)))
            new_edge_list = StatsBase.sample(params.edge_list, count_arcs_to_remove, replace=false)
            new_edge_count = length(new_edge_list)
        else
            new_edge_list = edges2remove
            new_edge_count = length(new_edge_list)
        end
        

        @info "Removed $(new_edge_count) our of $(length(params.edge_list)) edges"
        cost_of_expansion = base_cost * ones(new_edge_count)

        gparams.removed_edges = new_edge_list
        gparams.removed_edges_costs = cost_of_expansion

        if alg == "IFW"
            o = optimizer
        elseif alg == "NLMO-IFW" || alg == "NLMO-P"
            o = optimizer.sopt
        else
            @error "Invalid algorithm $alg"
        end

        add_network_design_constraints!(
            o,
            new_edge_list,
            con_list,
            g1, 
            params
        )

        if alg == "NLMO-IFW" || alg == "NLMO-P"
            for cidx in MOI.get(optimizer, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.ZeroOne}())
                push!(optimizer.integer_variables, cidx.value)
            end
        end
    end
    

    lmo = initialize_lmo(alg, optimizer, g1.graph, gparams)

    f, grad! = build_f_and_grad(alg, tn, gparams, cost_of_expansion, res)

    

    @info "Solving Problem!"

    od_pair_count = gparams.od_pair_count
    no_edges = length(gparams.init_nodes)

    #Choose correct algorithm depending on whether network design constraints are to be used or not. 
    if alg == "IFW"
        x, _, result = Boscia.solve(
            f, 
            grad!, 
            lmo, 
            verbose = bs_verb, 
            time_limit=time_limit, 
            max_fw_iter=max_fw_iter, 
            print_iter = 1,
            fw_verbose = fw_verb,
            dual_tightening=dual_tightening,
            rel_dual_gap=rel_dual_gap,
            line_search = FrankWolfe.Adaptive(verbose=ls_verb),
            use_postsolve = false
            #strong_convexity = minimum(cost_of_expansion)/2.001
            )

        obj = result[:primal_objective]
        non = result[:number_nodes]
        time_per_node = result[:total_time_in_sec] / non
        rel_dual_gap = result[:rel_dual_gap]

    elseif alg == "NLMO-FW"
        x0 = collect(FrankWolfe.compute_extreme_point(lmo, ones((od_pair_count + 1)* no_edges)))
        if any(x0 .< 0) == true && (maximum(x0[x0 .< 0]) < -1e-8)
            if maximum(x0[x0 .< 0]) < -1e-8
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x0[x0 .< 0] .= 0.0
        end

        x, _, obj = FrankWolfe.blended_pairwise_conditional_gradient(
                                            f, 
                                            grad!, 
                                            lmo, 
                                            copy(x0), 
                                            verbose=bs_verb, 
                                            lazy=true,
                                            line_search = FrankWolfe.Adaptive(verbose=ls_verb), 
                                            max_iteration=max_fw_iter,
                                            timeout=time_limit
                                            )
    elseif alg == "FW"
        x0 = collect(FrankWolfe.compute_extreme_point(lmo, ones((od_pair_count + 1)* no_edges)))
        if any(x0 .< 0) == true && (maximum(x0[x0 .< 0]) < -1e-8)
            if maximum(x0[x0 .< 0]) < -1e-8
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x0[x0 .< 0] .= 0.0
        end

        x, _, obj = FrankWolfe.blended_pairwise_conditional_gradient(
            f, 
            grad!, 
            lmo, 
            copy(x0), 
            verbose=bs_verb,
            lazy=true,
            line_search = FrankWolfe.Adaptive(verbose=ls_verb), 
            max_iteration=max_fw_iter,
            timeout=time_limit
            )
    elseif alg == "NLMO-IFW"
        optimizer.f = f

        x, _, result = Boscia.solve(
            f, 
            grad!, 
            lmo, 
            verbose=bs_verb,
            time_limit=time_limit, 
            max_fw_iter=max_fw_iter, 
            print_iter = 1, 
            dual_tightening=dual_tightening,
            fw_verbose = false,
            line_search=FrankWolfe.Adaptive(verbose=ls_verb),
            use_postsolve = false,
            #strong_convexity = minimum(cost_of_expansion)/2.001
            )

        #result = result[:primal_objective]
        obj = result[:primal_objective]
        non = result[:number_nodes]
        time_per_node = result[:total_time_in_sec] / non
        rel_dual_gap = result[:rel_dual_gap]
        
    else
        throw(ArgumentError("Incorrect Argument for method. $(alg)"))
    end

    if alg == "IFW" || alg == "NLMO-IFW"
        x_nd = x[((od_pair_count+1)*no_edges+1):((od_pair_count+1)*no_edges + new_edge_count)]
        new_edges_opened = sum(x_nd)
    end

    if alg == "IFW" || alg == "NLMO-IFW"
        #return [ta_data_name, alg, obj, base_cost, new_edges_opened, new_edge_count, x, new_edge_list, non]
        return Solution(ta_data_name, alg, obj, base_cost, new_edges_opened, 
                    new_edge_count, [], x, new_edge_list, 
                    non, time_per_node, rel_dual_gap)
    else
        #return [ta_data_name, alg, obj, base_cost, missing, missing]#, x, new_edge_list
        return Solution(ta_data_name, alg, obj, base_cost, 0, 
                    0, [], x, [], 
                    0, 0, 0)
    end
end