function all_or_nothing_single(travel_time_all::Vector{Float64}, params, graph)
    local state::Graphs.DijkstraState{Float64, Int}
    #x = zeros(size(params.init_nodes))

    edge_count = Graphs.ne(graph)
    od_pair_count = params.od_pair_count

    x = zeros(length(travel_time_all))
    edge_list = [(edge.src, edge.dst) for edge in collect(edges(graph))]
    @assert length(travel_time_all) == edge_count * od_pair_count + edge_count 

    travel_time = travel_time_all[edge_count* od_pair_count + 1:edge_count * (od_pair_count + 1)]     #take only the aggregated components
    
    init_nodes = [i for (i,j) in edge_list]
    term_nodes = [j for (i,j) in edge_list]
    new_link_dic = sparse(init_nodes, term_nodes, collect(1:edge_count))

    for r in 1:size(params.travel_demand)[1]
        # for each origin node r, find shortest paths to all destination nodes
        state = TA_dijkstra_shortest_paths(graph, travel_time, r, new_link_dic, params.first_thru_node)

        if any(iszero.(state.parents))
            unreachable_nodes = findall(iszero, state.dists)
            for ur in unreachable_nodes
                if ur != r
                    #println("Unreachable node: $(ur)")
                    if ur in 1:od_pair_count && params.travel_demand[r, ur] > 0.0
                        return MOI.INFEASIBLE
                    end
                end
            end
        end

        for s in 1:size(params.travel_demand)[2]
            # for each destination node s, find the shortest-path vector
            # load travel demand
            # x = x + travel_demand[r,s] * get_vector(state, r, s, link_dic)
            add_demand_vector_detailed!(x, params.travel_demand[r,s], state, r, s, new_link_dic, edge_list, od_pair_count)
        end

    end

    @assert all(x .>= 0)

    for edge in edge_list
        agg_x_flow = x[length(edge_list) * od_pair_count + new_link_dic[edge[1],edge[2]]]
        sum_x_flow = sum(x[new_link_dic[edge[1],edge[2]] + (j-1) * edge_count] for j in 1:od_pair_count)
        #println("edge: $(edge), agg_x_flow: $(agg_x_flow), sum_x_flow: $(sum_x_flow)")
        @assert agg_x_flow ≈ sum_x_flow
    end

    return x
end

function all_or_nothing(travel_time::Vector{Float64}, params, graph)
    # if nprocs() > 1 # if multiple CPU processes are available
    #     return all_or_nothing_parallel(travel_time, td, graph, link_dic)
    # else
    #     return all_or_nothing_single(travel_time, td, graph, link_dic)
    #     # when nprocs()==1, using @distributed just adds unnecessary setup time. I guess.
    # end

    return all_or_nothing_single(travel_time, params, graph)
end


function all_or_nothing_detailed(travel_time_all::Vector{Float64}, params, graph, true_graph)
    local state::Graphs.DijkstraState{Float64, Int}
    #x = zeros(size(params.init_nodes))

    edge_count = Graphs.ne(graph)
    od_pair_count = params.od_pair_count

    x = zeros(length(travel_time_all))
    edge_list = [(edge.src, edge.dst) for edge in collect(edges(graph))]

    println("length of travel_time_all = $(length(travel_time_all))")
    println("rhs = $( edge_count * od_pair_count + edge_count )")
    @assert length(travel_time_all) == edge_count * od_pair_count + edge_count 

    travel_time = travel_time_all[edge_count* od_pair_count + 1:edge_count * (od_pair_count + 1)]     #take only the aggregated components
    
    init_nodes = [i for (i,j) in edge_list]
    term_nodes = [j for (i,j) in edge_list]
    new_link_dic = sparse(init_nodes, term_nodes, collect(1:edge_count))

    dual_var_r = Dict([i=>zeros(od_pair_count) for i in range(1,step=1,stop=od_pair_count)])
    dual_var_t = zeros(od_pair_count)

    reverse_potentials = Dict()
    parents = Dict()

    for r in 1:size(params.travel_demand)[1]
        # for each origin node r, find shortest paths to all destination nodes
        state = TA_dijkstra_shortest_paths(graph, travel_time, r, new_link_dic, params.first_thru_node)

        if any(iszero.(state.parents))
            graph_unreachable_nodes = findall(iszero, state.parents)
            for ur in graph_unreachable_nodes
                if ur != r
                    if ur in 1:od_pair_count && params.travel_demand[r, ur] > 0.0
                        #println("Graph Unreachable node: $(ur)")
                        true_st = Graphs.dijkstra_shortest_paths(true_graph, r, allpaths=true)
                        return MOI.INFEASIBLE, state, r
                    end
                end
            end
        end

        if any(isinf.(state.dists))
            dist_unreachable_nodes = findall(isinf, state.dists)
            for ur in dist_unreachable_nodes
                if ur != r                    
                    if ur in 1:od_pair_count && params.travel_demand[r, ur] > 0.0 && isinf(state.dists[state.parents[ur]])
                        #println("Distance Unreachable node: $(ur)")
                        return MOI.INFEASIBLE, state, r
                    end
                end
            end
        end

        reverse_potentials[r] = state.dists
        parents[r] = state.parents
        
        for s in 1:size(params.travel_demand)[2]
            # for each destination node s, find the shortest-path vector
            # load travel demand
            # x = x + travel_demand[r,s] * get_vector(state, r, s, link_dic)

            #the case of infinite cost means that the path travels through a zonal node which is invalid. 
            if state.parents[s] != 0 && (isinf(state.dists[state.parents[s]]) == false)
                #@assert isinf(state.dists[state.parents[s]]) == false
                dual_var_r[r][s] = state.dists[state.parents[s]] + travel_time[new_link_dic[state.parents[s] ,s]] #distance from source r to destination s
            end

            add_demand_vector_detailed!(x, params.travel_demand[r,s], state, r, s, new_link_dic, edge_list, od_pair_count)
        end

    end

    @assert all(x .>= 0)

    for edge in edge_list
        agg_x_flow = x[length(edge_list) * od_pair_count + new_link_dic[edge[1],edge[2]]]
        sum_x_flow = sum(x[new_link_dic[edge[1],edge[2]] + (j-1) * edge_count] for j in 1:od_pair_count)
        #println("edge: $(edge), agg_x_flow: $(agg_x_flow), sum_x_flow: $(sum_x_flow)")
        @assert agg_x_flow ≈ sum_x_flow
    end

    return x, dual_var_r, dual_var_t, reverse_potentials, parents
end


function all_or_nothing_LP(scen, travel_time_all::Vector{Float64}, params, y_val, node_count, edge_list, GRB_ENV)

    od_pair_count = params.od_pair_count


    #edge_list = [(params.init_nodes[i], params.term_nodes[i]) for i in 1:length(params.init_nodes)]
    edge_count = length(edge_list)

    travel_time = travel_time_all[edge_count* od_pair_count + 1:edge_count * (od_pair_count + 1)]

    first_non_zone_node = params.first_thru_node

    #demands = params.travel_demand[scen]

    

    arc_cost = 1e8 * ones(node_count, node_count)

    for (idx, (i,j)) in enumerate(edge_list)
        arc_cost[i,j] = travel_time[params.link_dic[i,j]]

        # if j < first_non_zone_node
        #     arc_cost[i,j] = 1e6
        # end
    end

    zone_outgoing_edges = [(i,j) for (i,j) in edge_list if i < od_pair_count + 1]

    #@info zone_outgoing_edges

    zone_incoming_edges = [(i,j) for (i,j) in edge_list if j < od_pair_count + 1]

    #@info zone_incoming_edges
    
    non_zone_edges = [(i,j) for (i,j) in edge_list if i >= od_pair_count + 1 && j >= od_pair_count + 1]

    removed_edges = params.removed_edges
    removed_edges_count = length(removed_edges)


    function get_edge_index(e)
        return findfirst(isequal(e), removed_edges)
    end

    nr_zone_outgoing_edges = [e for e in zone_outgoing_edges if !(e in removed_edges)]
    nr_zone_incoming_edges = [e for e in zone_incoming_edges if !(e in removed_edges)]
    nr_non_zone_edges = [e for e in non_zone_edges if !(e in removed_edges)]

    r_zone_outgoing_edges = [e for e in zone_outgoing_edges if (e in removed_edges)]
    r_zone_incoming_edges = [e for e in zone_incoming_edges if (e in removed_edges)]
    r_non_zone_edges = [e for e in non_zone_edges if (e in removed_edges)]

    model = Model(() -> Gurobi.Optimizer(GRB_ENV)) 
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "InfUnbdInfo", 1)

    @variable(model, q[1:node_count, 1:od_pair_count])
    @variable(model, s[1:removed_edges_count, 1:od_pair_count] >= 0)
    @variable(model, r[1:od_pair_count, 1:od_pair_count])
    @variable(model, t[1:od_pair_count]>=0)

    #@info "Added $(node_count * od_pair_count) variables (q)"
    #@info "Added $(removed_edges_count * od_pair_count) variables (s)"
    #@info "Added $(od_pair_count * od_pair_count) variables (r)"
    #@info "Added $(od_pair_count) variables (t)"

    @constraint(
        model, 
        potential_across_non_zone_arcs[(i,j) in nr_non_zone_edges, z in 1:od_pair_count],
        q[i,z] - q[j,z] <= arc_cost[i,j]
        )
    
    #@info "Added $(length(nr_non_zone_edges)*od_pair_count) potential constraints (non-zone arcs)"

    @constraint(
        model, 
        potential_across_outgoing_arcs[(i,j) in nr_zone_outgoing_edges, z in 1:od_pair_count],
        r[i,z] - q[j,z] <= arc_cost[i,j]
        )

    #@info "Added $(length(nr_zone_outgoing_edges)*od_pair_count) potential constraints (outgoing arcs)"

    @constraint(
        model, 
        potential_across_incoming_arcs[(i,z) in nr_zone_incoming_edges],
        q[i,z] - t[z] <= arc_cost[i,z]
        )

    #@info "Added $(length(nr_zone_incoming_edges)) potential constraints (incoming arcs)"
        
    @constraint(
        model, 
        potential_across_removed_non_zone_arcs[(i,j) in r_non_zone_edges, z in 1:od_pair_count],
        q[i,z] - q[j,z] - s[get_edge_index((i,j)),z] <= arc_cost[i,j]
        )

    #@info "Added $(length(r_non_zone_edges)*od_pair_count) potential constraints (removed non zone arcs)"

    @constraint(
        model, 
        potential_across_removed_outgoing_arcs[(i,j) in r_zone_outgoing_edges, z in 1:od_pair_count],
        r[i,z] - q[j,z] - s[get_edge_index((i,j)),z] <= arc_cost[i,j]
        )
        
    #@info "Added $(length(r_zone_outgoing_edges)*od_pair_count) potential constraints (removed outgoing arcs)"
        
    @constraint(
        model, 
        potential_across_removed_incoming_arcs[(i,z) in r_zone_incoming_edges],
        q[i,z] - t[z] - s[get_edge_index((i,z)),z] <= arc_cost[i,z]
        )

    #@info "Added $(length(r_zone_incoming_edges)) potential constraints (removed incoming arcs)"

    tot_constraints = length(r_zone_incoming_edges) + length(r_zone_outgoing_edges)*od_pair_count + length(r_non_zone_edges)*od_pair_count 
                    + length(nr_zone_incoming_edges) + length(nr_zone_outgoing_edges)*od_pair_count + length(nr_non_zone_edges)*od_pair_count

    #@info "Added total of $(tot_constraints) potential constraints"

    r_non_zone_edges

    td_2_dest = [maximum([sum(params.travel_demand[s][:,z]) for s in 1:params.total_scenarios]) for z in 1:params.od_pair_count]

    #td_2_dest = [sum(demands[:, z]) for z in 1:od_pair_count] #Total travel demand to each destination

    travel_cost_expression = @expression(
        model,
        sum(
            sum((r[i,z] - t[z]) * params.travel_demand[scen][i,z] for i in 1:od_pair_count) 
            for z in 1:od_pair_count
            )  
        )

    removed_edge_expression = @expression(
        model,
        sum(
            td_2_dest[z] * sum(s[get_edge_index((i,j)),z] * y_val[get_edge_index((i,j))] for (i,j) in removed_edges) 
            for z in 1:od_pair_count
            )
        )

    @objective(
        model,
        Max,
        travel_cost_expression - removed_edge_expression
        )

    #write_to_file(model, "model5.lp")
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        dual_var_r = value.(r)
        dual_var_t = value.(t)
        dual_var_s = value.(s)
        return MOI.INFEASIBLE, dual_var_r, dual_var_t, dual_var_s 
    end

    dual_var_r = value.(r)
    dual_var_t = value.(t)
    dual_var_s = value.(s)

    return MOI.OPTIMAL, dual_var_r, dual_var_t, dual_var_s, value.(q)
end


function all_or_nothing_LP_no_zones(travel_time_all::Vector{Float64}, params, y_val, node_count, GRB_ENV)

    od_pair_count = params.od_pair_count


    edge_list = [(params.init_nodes[i], params.term_nodes[i]) for i in 1:length(params.init_nodes)]
    edge_count = length(edge_list)

    travel_time = travel_time_all[edge_count*od_pair_count+1:edge_count*(od_pair_count+1)]

    first_non_zone_node = params.first_thru_node

    demands = params.travel_demand



    arc_cost = 1e6 * ones(node_count, node_count)

    for (idx, (i, j)) in enumerate(edge_list)
        arc_cost[i, j] = travel_time[params.link_dic[i, j]]

        # if j < first_non_zone_node
        #     arc_cost[i,j] = 1e6
        # end
    end

    zone_outgoing_edges = [(i, j) for (i, j) in edge_list if i < od_pair_count + 1]

    #@info zone_outgoing_edges

    zone_incoming_edges = [(i, j) for (i, j) in edge_list if j < od_pair_count + 1]

    #@info zone_incoming_edges

    non_zone_edges = [(i, j) for (i, j) in edge_list if i >= od_pair_count + 1 && j >= od_pair_count + 1]

    removed_edges = params.removed_edges
    removed_edges_count = length(removed_edges)


    function get_edge_index(e)
        return findfirst(isequal(e), removed_edges)
    end

    nr_zone_outgoing_edges = [e for e in zone_outgoing_edges if !(e in removed_edges)]
    nr_zone_incoming_edges = [e for e in zone_incoming_edges if !(e in removed_edges)]
    nr_non_zone_edges = [e for e in non_zone_edges if !(e in removed_edges)]

    r_zone_outgoing_edges = [e for e in zone_outgoing_edges if (e in removed_edges)]
    r_zone_incoming_edges = [e for e in zone_incoming_edges if (e in removed_edges)]
    r_non_zone_edges = [e for e in non_zone_edges if (e in removed_edges)]

    model = Model(() -> Gurobi.Optimizer(GRB_ENV))
    set_optimizer_attribute(model, "OutputFlag", 0)
    set_optimizer_attribute(model, "InfUnbdInfo", 1)

    @variable(model, q[1:node_count, 1:od_pair_count])
    @variable(model, s[1:removed_edges_count, 1:od_pair_count] >= 0)
    @variable(model, r[1:od_pair_count, 1:od_pair_count])
    @variable(model, t[1:od_pair_count])

    #@info "Added $(node_count * od_pair_count) variables (q)"
    #@info "Added $(removed_edges_count * od_pair_count) variables (s)"
    #@info "Added $(od_pair_count * od_pair_count) variables (r)"
    #@info "Added $(od_pair_count) variables (t)"

    @constraint(
        model,
        potential_across_arcs[(i, j) in edge_list, z in 1:od_pair_count],
        q[i, z] - q[j, z] <= arc_cost[i, j]
    )

    #@info "Added $(length(nr_non_zone_edges)*od_pair_count) potential constraints (non-zone arcs)"

    @constraint(
        model,
        potential_across_outgoing_arcs[(i, j) in nr_zone_outgoing_edges, z in 1:od_pair_count],
        r[i, z] - q[j, z] <= arc_cost[i, j]
    )

    #@info "Added $(length(nr_zone_outgoing_edges)*od_pair_count) potential constraints (outgoing arcs)"

    @constraint(
        model,
        potential_across_incoming_arcs[(i, z) in nr_zone_incoming_edges],
        q[i, z] - t[z] <= arc_cost[i, z]
    )

    #@info "Added $(length(nr_zone_incoming_edges)) potential constraints (incoming arcs)"

    @constraint(
        model,
        potential_across_removed_non_zone_arcs[(i, j) in r_non_zone_edges, z in 1:od_pair_count],
        q[i, z] - q[j, z] - s[get_edge_index((i, j)), z] <= arc_cost[i, j]
    )

    #@info "Added $(length(r_non_zone_edges)*od_pair_count) potential constraints (removed non zone arcs)"

    @constraint(
        model,
        potential_across_removed_outgoing_arcs[(i, j) in r_zone_outgoing_edges, z in 1:od_pair_count],
        r[i, z] - q[j, z] - s[get_edge_index((i, j)), z] <= arc_cost[i, j]
    )

    #@info "Added $(length(r_zone_outgoing_edges)*od_pair_count) potential constraints (removed outgoing arcs)"

    @constraint(
        model,
        potential_across_removed_incoming_arcs[(i, z) in r_zone_incoming_edges],
        q[i, z] - t[z] - s[get_edge_index((i, z)), z] <= arc_cost[i, z]
    )

    #@info "Added $(length(r_zone_incoming_edges)) potential constraints (removed incoming arcs)"

    tot_constraints = length(r_zone_incoming_edges) + length(r_zone_outgoing_edges) * od_pair_count + length(r_non_zone_edges) * od_pair_count
    length(nr_zone_incoming_edges) + length(nr_zone_outgoing_edges) * od_pair_count + length(nr_non_zone_edges) * od_pair_count

    #@info "Added total of $(tot_constraints) potential constraints"

    +

    r_non_zone_edges

    td_2_dest = [sum(demands[:, z]) for z in 1:od_pair_count] #Total travel demand to each destination

    travel_cost_expression = @expression(
        model,
        sum(
            sum((r[i, z] - t[z]) * demands[i, z] for i in 1:od_pair_count)
            for z in 1:od_pair_count
        )
    )

    removed_edge_expression = @expression(
        model,
        sum(
            td_2_dest[z] * sum(s[get_edge_index((i, j)), z] * y_val[get_edge_index((i, j))] for (i, j) in removed_edges)
            for z in 1:od_pair_count
        )
    )

    @objective(
        model,
        Max,
        travel_cost_expression - removed_edge_expression
    )

    #write_to_file(model, "model5.lp")
    optimize!(model)

    if termination_status(model) != MOI.OPTIMAL
        dual_var_r = value.(r)
        dual_var_t = value.(t)
        dual_var_s = value.(s)
        return MOI.INFEASIBLE, dual_var_r, dual_var_t, dual_var_s
    end

    dual_var_r = value.(r)
    dual_var_t = value.(t)
    dual_var_s = value.(s)

    return MOI.OPTIMAL, dual_var_r, dual_var_t, dual_var_s
end


function all_or_nothing_stochastic_detailed(travel_time_all::Vector{Float64}, params, graph, true_graph, scen, distmx)
    local state::Graphs.DijkstraState{Float64,Int}
    #x = zeros(size(params.init_nodes))

    edge_count = Graphs.ne(graph)
    od_pair_count = params.od_pair_count

    x = zeros(length(travel_time_all))
    edge_list = [(edge.src, edge.dst) for edge in collect(edges(graph))]
    #println("length of travel_time_all = $(length(travel_time_all))")
    #println("rhs = $( edge_count * od_pair_count + edge_count )")
    @assert length(travel_time_all) == edge_count * od_pair_count + edge_count

    travel_time = travel_time_all[edge_count*od_pair_count+1:edge_count*(od_pair_count+1)]     #take only the aggregated components

    init_nodes = [i for (i, j) in edge_list]
    term_nodes = [j for (i, j) in edge_list]
    new_link_dic = sparse(init_nodes, term_nodes, collect(1:edge_count))

    dual_var_r = Dict([i => zeros(od_pair_count) for i in range(1, step=1, stop=od_pair_count)])
    dual_var_t = zeros(od_pair_count)

    reverse_potentials = Dict()
    parents = Dict()

    for r in 1:size(params.travel_demand[scen])[1]
        # for each origin node r, find shortest paths to all destination nodes
        state = TA_dijkstra_shortest_paths(graph, travel_time, r, new_link_dic, params.first_thru_node, distmx)

        if any(iszero.(state.parents))
            graph_unreachable_nodes = findall(iszero, state.parents)
            for ur in graph_unreachable_nodes
                if ur != r
                    if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0
                        #println("Graph Unreachable node: $(ur)")
                        true_st = Graphs.dijkstra_shortest_paths(true_graph, r, allpaths=true)
                        return MOI.INFEASIBLE, state, r
                    end
                end
            end
        end

        if any(isinf.(state.dists))
            dist_unreachable_nodes = findall(isinf, state.dists)
            for ur in dist_unreachable_nodes
                if ur != r
                    if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0 && isinf(state.dists[state.parents[ur]])
                        #println("Distance Unreachable node: $(ur)")
                        return MOI.INFEASIBLE, state, r
                    end
                end
            end
        end

        reverse_potentials[r] = state.dists
        parents[r] = state.parents

        for s in 1:size(params.travel_demand[scen])[2]
            # for each destination node s, find the shortest-path vector
            # load travel demand
            # x = x + travel_demand[r,s] * get_vector(state, r, s, link_dic)

            #the case of infinite cost means that the path travels through a zonal node which is invalid. 
            if state.parents[s] != 0 && (isinf(state.dists[state.parents[s]]) == false)
                #@assert isinf(state.dists[state.parents[s]]) == false
                #println("r: $(r), s: $(s)")
                #println("state.parents[s]: $(state.parents[s])")
                #println("new_link_dic[state.parents[s], s]: $(new_link_dic[state.parents[s], s])")
                #println("travel_time[new_link_dic[state.parents[s], s]]: $(travel_time[new_link_dic[state.parents[s], s]])")
                dual_var_r[r][s] = state.dists[state.parents[s]] + travel_time[new_link_dic[state.parents[s], s]] #distance from source r to destination s
            end

            #add_demand_vector_detailed!(x, params.travel_demand[r, s], state, r, s, new_link_dic, edge_list, od_pair_count)
            add_demand_vector_detailed!(x, params.travel_demand[scen][r, s], state, r, s, new_link_dic, edge_list, od_pair_count)
        end

    end

    @assert all(x .>= 0)

    for edge in edge_list
        agg_x_flow = x[length(edge_list)*od_pair_count+new_link_dic[edge[1], edge[2]]]
        sum_x_flow = sum(x[new_link_dic[edge[1], edge[2]]+(j-1)*edge_count] for j in 1:od_pair_count)
        #println("edge: $(edge), agg_x_flow: $(agg_x_flow), sum_x_flow: $(sum_x_flow)")
        @assert agg_x_flow ≈ sum_x_flow
    end

    return x, dual_var_r, dual_var_t, reverse_potentials, parents
end


function all_or_nothing_stochastic_separated(travel_time_all::Vector{Float64}, params, graph, curr_edge_list, scen, distmx)
    local state::Graphs.DijkstraState{Float64,Int}
    #x = zeros(size(params.init_nodes))

    edge_count = Graphs.ne(graph)
    od_pair_count = params.od_pair_count

    x = zeros(length(travel_time_all))
    #edge_list = [(edge.src, edge.dst) for edge in collect(edges(graph))]
    @assert length(travel_time_all) == edge_count * od_pair_count + edge_count

    init_nodes = [i for (i, j) in curr_edge_list]
    term_nodes = [j for (i, j) in curr_edge_list]
    new_link_dic = sparse(init_nodes, term_nodes, collect(1:edge_count))

    for r in 1:size(params.travel_demand[scen])[1]
        for s in 1:size(params.travel_demand[scen])[2]

            travel_time = travel_time_all[(s-1) * edge_count+1:s * edge_count]     #take only the aggregated components

            # for each origin node r and s pair, find shortest paths to all destination nodes
            state = TA_dijkstra_shortest_paths(graph, travel_time, r, new_link_dic, params.first_thru_node, distmx)
            # for each destination node s, find the shortest-path vector
            # load travel demand
            # x = x + travel_demand[r,s] * get_vector(state, r, s, link_dic)


            if any(iszero.(state.parents))
                graph_unreachable_nodes = findall(iszero, state.parents)
                for ur in graph_unreachable_nodes
                    if ur != r
                        if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0
                            return MOI.INFEASIBLE
                        end
                    end
                end
            end

            if any(isinf.(state.dists))
                dist_unreachable_nodes = findall(isinf, state.dists)
                for ur in dist_unreachable_nodes
                    if ur != r
                        if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0 && isinf(state.dists[state.parents[ur]])
                            return MOI.INFEASIBLE
                        end
                    end
                end
            end

            #add_demand_vector_detailed!(x, params.travel_demand[r, s], state, r, s, new_link_dic, edge_list, od_pair_count)
            add_demand_vector_detailed!(x, params.travel_demand[scen][r, s], state, r, s, new_link_dic, curr_edge_list, od_pair_count)
        end

    end

    @assert all(x .>= 0)

    for edge in curr_edge_list
        agg_x_flow = x[length(curr_edge_list)*od_pair_count+new_link_dic[edge[1], edge[2]]]
        sum_x_flow = sum(x[new_link_dic[edge[1], edge[2]]+(j-1)*edge_count] for j in 1:od_pair_count)
        #println("edge: $(edge), agg_x_flow: $(agg_x_flow), sum_x_flow: $(sum_x_flow)")
        @assert agg_x_flow ≈ sum_x_flow
    end

    return x
end

#costs are duplicated for all destinations. 
function all_or_nothing_stochastic_common(travel_time_all::Vector{Float64}, params, graph, true_graph, scen, distmx)
    local state::Graphs.DijkstraState{Float64,Int}
    #x = zeros(size(params.init_nodes))

    edge_count = Graphs.ne(graph)
    od_pair_count = params.od_pair_count

    x = zeros(length(travel_time_all))
    edge_list = [(edge.src, edge.dst) for edge in collect(edges(graph))]
    @assert length(travel_time_all) == edge_count * od_pair_count + edge_count

    init_nodes = [i for (i, j) in edge_list]
    term_nodes = [j for (i, j) in edge_list]
    new_link_dic = sparse(init_nodes, term_nodes, collect(1:edge_count))

    travel_time = travel_time_all[1:edge_count]     

    for r in 1:size(params.travel_demand[scen])[1]
        # for each origin node r and s pair, find shortest paths to all destination nodes
        state = TA_dijkstra_shortest_paths(graph, travel_time, r, new_link_dic, params.first_thru_node, distmx)

        if any(iszero.(state.parents))
            graph_unreachable_nodes = findall(iszero, state.parents)
            for ur in graph_unreachable_nodes
                if ur != r
                    if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0
                        return MOI.INFEASIBLE
                    end
                end
            end
        end

        if any(isinf.(state.dists))
            dist_unreachable_nodes = findall(isinf, state.dists)
            for ur in dist_unreachable_nodes
                if ur != r
                    if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0 && isinf(state.dists[state.parents[ur]])
                        return MOI.INFEASIBLE
                    end
                end
            end
        end

        for s in 1:size(params.travel_demand[scen])[2]

            
            # for each destination node s, find the shortest-path vector
            # load travel demand
            # x = x + travel_demand[r,s] * get_vector(state, r, s, link_dic)


            #add_demand_vector_detailed!(x, params.travel_demand[r, s], state, r, s, new_link_dic, edge_list, od_pair_count)
            add_demand_vector_detailed!(x, params.travel_demand[scen][r, s], state, r, s, new_link_dic, edge_list, od_pair_count)
        end

    end

    @assert all(x .>= 0)

    for edge in edge_list
        agg_x_flow = x[length(edge_list)*od_pair_count+new_link_dic[edge[1], edge[2]]]
        sum_x_flow = sum(x[new_link_dic[edge[1], edge[2]]+(j-1)*edge_count] for j in 1:od_pair_count)
        #println("edge: $(edge), agg_x_flow: $(agg_x_flow), sum_x_flow: $(sum_x_flow)")
        @assert agg_x_flow ≈ sum_x_flow
    end

    return x
end

function all_or_nothing_stochastic_reversed(travel_time_all::Vector{Float64}, params, graph, curr_edge_list, scen, distmx)
    local state::Graphs.DijkstraState{Float64,Int}
    #x = zeros(size(params.init_nodes))

    edge_count = Graphs.ne(graph)
    od_pair_count = params.od_pair_count

    x = zeros(length(travel_time_all))
    #edge_list = [(edge.src, edge.dst) for edge in collect(edges(graph))]
    @assert length(travel_time_all) == edge_count * od_pair_count + edge_count

    init_nodes = [i for (i, j) in curr_edge_list]
    term_nodes = [j for (i, j) in curr_edge_list]


    new_link_dic = sparse(term_nodes, init_nodes, collect(1:edge_count))

    for s in 1:size(params.travel_demand[scen])[2]

        if sum(params.travel_demand[scen][:,s]) == 0.0
            continue
        end

        travel_time = travel_time_all[(s-1) * edge_count+1:s * edge_count]     
        
        # for each origin node r and s pair, find shortest paths to all destination nodes
        
        state = TA_reverse_shortest_paths(graph, travel_time, s, new_link_dic, params.first_thru_node, distmx)

        
        for r in 1:size(params.travel_demand[scen])[1]

            
            # for each destination node s, find the shortest-path vector
            # load travel demand
            # x = x + travel_demand[r,s] * get_vector(state, r, s, link_dic)

            if state.parents[r] == 0 && params.travel_demand[scen][r, s] > 0.0
                return MOI.INFEASIBLE
            end

            # if any(iszero.(state.parents))
            #     graph_unreachable_nodes = findall(iszero, state.parents)
            #     for ur in graph_unreachable_nodes
            #         if ur != r
            #             if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0
            #                 println("r: $(r), s: $(s), ur: $(ur)")
            #                 println("travel demand: $(params.travel_demand[scen][r, ur])")
            #                 return MOI.INFEASIBLE
            #             end
            #         end
            #     end
            # end
            

            # if any(isinf.(state.dists))
            #     dist_unreachable_nodes = findall(isinf, state.dists)
            #     for ur in dist_unreachable_nodes
            #         if ur != r
            #             if ur in 1:od_pair_count && params.travel_demand[scen][r, ur] > 0.0 && isinf(state.dists[state.parents[ur]])
            #                 return MOI.INFEASIBLE
            #             end
            #         end
            #     end
            # end

            #for the reverse case we have flipped the source and the destination
            #add_demand_vector_detailed!(x, params.travel_demand[r, s], state, r, s, new_link_dic, edge_list, od_pair_count)
            add_demand_vector_detailed_reverse!(x, params.travel_demand[scen][r, s], state, s, r, new_link_dic, curr_edge_list, od_pair_count)
        end

    end

    @assert all(x .>= 0)

    for edge in curr_edge_list
        agg_x_flow = x[length(curr_edge_list)*od_pair_count+new_link_dic[edge[2], edge[1]]]
        sum_x_flow = sum(x[new_link_dic[edge[2], edge[1]]+(j-1)*edge_count] for j in 1:od_pair_count)
        #println("edge: $(edge), agg_x_flow: $(agg_x_flow), sum_x_flow: $(sum_x_flow)")
        @assert agg_x_flow ≈ sum_x_flow
    end

    return x
end

function create_reverse_graph(graph)
    reverse_graph = Graphs.SimpleDiGraph(Graphs.nv(graph))
    for edge in collect(edges(graph))
        Graphs.add_edge!(reverse_graph, edge.dst, edge.src)
    end
    return reverse_graph
end