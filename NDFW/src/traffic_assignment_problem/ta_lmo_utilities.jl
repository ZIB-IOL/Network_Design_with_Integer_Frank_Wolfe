function initialize_lmo_ta(alg, optimizer, g, gparams)
    if alg == "LP"
        lmo = FrankWolfe.MathOptLMO(optimizer, false)
    elseif alg == "SP"
        lmo = SPLMO(g, gparams)
    else
        throw(ArgumentError("Incorrect Argument for method. $(alg)"))
    end

    return lmo
end

#costs are duplicated for all destinations. 
function all_or_nothing_sp(travel_time::Vector{Float64}, params, graph, scen, link_dic, distmx)
    local state::Graphs.DijkstraState{Float64,Int}
    #x = zeros(size(params.init_nodes))


    x = zeros(length(travel_time))

    for r in 1:size(params.travel_demand[scen])[1]
        # for each origin node r and s pair, find shortest paths to all destination nodes
        #state = TA_dijkstra_shortest_paths(graph, travel_time, r, link_dic, params.first_thru_node, distmx)
        state = TA_dijkstra_shortest_paths_default(graph, travel_time, r, params.init_nodes, params.term_nodes, params.first_thru_node, distmx)
        
        for s in 1:size(params.travel_demand[scen])[2]


            # for each destination node s, find the shortest-path vector
            # load travel demand
            # x = x + travel_demand[r,s] * get_vector(state, r, s, link_dic)


            #add_demand_vector_detailed!(x, params.travel_demand[r, s], state, r, s, new_link_dic, edge_list, od_pair_count)
            add_demand_vector!(x, params.travel_demand[scen][r, s], state, r, s, link_dic)
        end

    end

    # @assert all(x .>= 0)

    # for edge in edge_list
    #     agg_x_flow = x[length(edge_list)*od_pair_count+link_dic[edge[1], edge[2]]]
    #     sum_x_flow = sum(x[link_dic[edge[1], edge[2]]+(j-1)*edge_count] for j in 1:od_pair_count)
    #     #println("edge: $(edge), agg_x_flow: $(agg_x_flow), sum_x_flow: $(sum_x_flow)")
    #     @assert agg_x_flow ≈ sum_x_flow
    # end

    return x
end