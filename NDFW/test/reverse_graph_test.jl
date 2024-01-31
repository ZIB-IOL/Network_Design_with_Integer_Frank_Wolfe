using NDFW
using Graphs
using SparseArrays

# edge_list = [(1,5), (2, 6), (5,6), (5,4), (6,3), (3,5)]
# ftn = 5

# graph = Graphs.SimpleDiGraph{Int64}(6)

# for edge in edge_list
#     Graphs.add_edge!(graph, edge[1], edge[2])
# end

# init_nodes = []
# term_nodes = []

# for edge in edge_list
#     push!(init_nodes, edge[1])
#     push!(term_nodes, edge[2])
# end

# costs = [1, 1, 1, 0, 0, 1]



# link_dic = sparse(init_nodes, term_nodes, collect(1:Graphs.ne(graph)))

# rgraph = NDFW.create_reverse_graph(graph)

# origin = 1
# dest = 2
# stf1 = NDFW.TA_dijkstra_shortest_paths(graph, costs, 1, link_dic, ftn)
# stf2 = NDFW.TA_dijkstra_shortest_paths(graph, costs, 2, link_dic, ftn)
# stb3 = NDFW.TA_reverse_shortest_paths(rgraph, costs, 3, init_nodes, term_nodes, ftn)
# stb4 = NDFW.TA_reverse_shortest_paths(rgraph, costs, 4, init_nodes, term_nodes, ftn)

ta_data_name = "Berlin-Tiergarten"
alg = "BNDLMO-PR"
base_cost = 1.0
remove_circular_demand = true
total_scenarios = 1
scale = 0.1

tn, g1, params, gparams = NDFW.initialize_parameters(ta_data_name, 
                            alg, 
                            base_cost, 
                            remove_circular_demand, 
                            total_scenarios, 
                            scale)

generated_points = false

initpt = nothing
destpt = nothing

# while generated_points == false
#     initpt = rand(1:tn.first_thru_node-1)
#     destpt = rand(1:tn.first_thru_node-1)
#     if initpt == destpt
#         continue
#     end
#     if gparams.travel_demand[1][initpt, destpt] == 0
#         continue
#     end
#     generated_points = true
# end

initpt = 1
destpt = 4

x = output.soln_vector
storage = zeros(var_count)
gblmo(storage, x)
gb = copy(storage)

#costs = rand(length(tn.init_node))
rg1 = NDFW.create_reverse_graph(g1.graph)

rev_link_dic = sparse(tn.term_node, tn.init_node, collect(1:Graphs.ne(rg1)))

distmx = zeros(Graphs.nv(g1.graph), Graphs.nv(g1.graph))

for destpt in 1:gparams.od_pair_count
    costs = gb[(destpt-1)*edge_count+1:destpt*edge_count]
    for initpt in 1:gparams.od_pair_count
        println("initpt:$(initpt), dsetpt:$(destpt)")

        stdfinit = NDFW.TA_dijkstra_shortest_paths(g1.graph, costs, initpt, gparams.link_dic, tn.first_thru_node, distmx)
        #stdbdest = NDFW.TA_dijkstra_shortest_paths(rg1, costs, destpt, rev_link_dic, tn.first_thru_node)

        stdbdest = NDFW.TA_reverse_shortest_paths(rg1, costs, destpt, rev_link_dic, tn.first_thru_node, distmx)

        pathid = NDFW.get_path(stdfinit, initpt, destpt)
        pathdi = NDFW.get_path(stdbdest, destpt, initpt)

        if pathid[1] != [(j, i) for (i, j) in reverse(pathdi[1])]
            println(pathid[1])
            println([(j, i) for (i, j) in reverse(pathdi[1])])
        end

    end
end