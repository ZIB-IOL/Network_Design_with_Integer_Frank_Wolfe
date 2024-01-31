using JuMP
using DataFrames
using NDFW
using Random
using Graphs
using Gurobi
using StatsBase
using Graphs
using SparseArrays


"""
This function use the dual variables for the source constraint and the distance from the source
to each node to calculate the potentials. 
"""
function get_potentials_fw_mod(g, travel_time, edge_list, ftn)

    init_nodes = [i for (i,j) in edge_list]
    term_nodes = [j for (i,j) in edge_list]
    new_link_dic = sparse(init_nodes, term_nodes, collect(1:length(edge_list)))


    state = NDFW.TA_floyd_warshall_shortest_paths(g, travel_time, new_link_dic, ftn)
    return state.dists, state.parents
end
    
"""
The function computes the distance from every origin node to a destination node. 
"""
function compute_zonal_distances(dists, parents, fft, edge_list, od_pair_count)
    zonal_distances = zeros(od_pair_count, od_pair_count)
    for i in 1:od_pair_count
        for j in 1:od_pair_count
            if (i != j)
                #@info "i = $i, j = $j"
                last_edge = findall((x)->x==(parents[i,j],j), edge_list)
                if length(last_edge) > 0
                    final_edge_idx = last_edge[1]
                    zonal_distances[i,j] = dists[i,parents[i,j]] + fft[final_edge_idx]
                else
                    zonal_distances[i,j] = 1e6
                end
            end
        end
    end
    return zonal_distances
end

function compute_potentials(dists, parents, node_count, 
                            od_pair_count, fft, edge_list, arc_cost, 
                            in_neighbours, out_neighbours)

    potentials = zeros(node_count, node_count)

    for i in 1:node_count
        for z in 1:od_pair_count
            #Destination node. By default the potential should be 0. 
            #We should probably verify this. 
            if i == z
                continue
            end

            #Current node
            @info "i = $i, z = $z"

            #By the construction of the graph, the distance to entering the graph is at least 1e6 because
            #this is the cost of the final arc. 
            if dists[i,z] < 2e6

                #index of parent arc. parents[i,z] is the parent of z in the path from i to z
                final_edge_idx = findall(x->x==(parents[i,z],z), edge_list)[1]

                #Distance from current node i to destination z is equal to the distance returned by 
                #the algorithm - the modified cost of final arc + actual cost of final arc.

                potentials[i] = dists[i,z] - 1e6 + fft[final_edge_idx]
            elseif dists[i,z] > 2e6
                # If it is not possible to reach the destination node from the current node, then the
                # potential is equal to 
                #       minimum (maximum of incoming potentials, minimum of outgoing potentials) - 1e6

                dinz = 0
                dotz = 0

                if length(in_neighbours[i]) > 0
                    all_inc = []
                    for inc in in_neighbours[i]
                        final_edge_idx = findall(x->x==(parents[inc,z],z), edge_list)[1]
                        diz = dists[inc,parents[inc,z]] + fft[final_edge_idx] - arc_cost[inc,i]
                        push!(all_inc, diz)
                    end

                    dinz = maximum(all_inc)
                end

                if length(out_neighbours[i]) > 0

                    all_out = []
                    for out in out_neighbours[i]
                        final_edge_idx = findall(x->x==(parents[out,z],z), edge_list)[1]
                        dotz = dists[out, parents[out, z]] + fft[final_edge_idx] + arc_cost[i, out]
                        push!(all_out, dotz)
                    end

                    dotz = minimum(all_out)
                end
                
                #I removed this since out going is not going to lead to anywhere 
                #potentials[i] = minimum(dinz, dotz) - 1e6
                @assert dinz <= dotz

                potentials[i] = dinz - 1e6

            end
        end
    end

    return potentials
end


function compute_s_vals(potentials, removed_edges, od_pair_count, arc_cost)
    s_vals = zeros(length(removed_edges), od_pair_count)

    for (idx, (i,j)) in enumerate(removed_edges)
        for z in 1:od_pair_count
            s_vals[idx, z] = max(potentials[i,z] - potentials[j,z] - arc_cost[i,j], 0)
        end
    end
    return s_vals
end


seed = 0
@show seed
Random.seed!(seed)

#ds = "Berlin-Mitte-Center"
ds = "Braess-Example"

tn = NDFW.load_ta_network(ds)

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

od_pair_count = zone_count > 1 ? zone_count : node_count
first_non_zone_node = tn.first_thru_node

demands = tn.travel_demand

edge_list = [(tn.init_node[i], tn.term_node[i]) for i in 1:edge_count]

g = Graphs.SimpleDiGraph(node_count)
for edge in edge_list
    add_edge!(g, edge[1], edge[2])
end


add_edge!(g, 3, 5)
add_edge!(g, 1, 5)
add_edge!(g, 5, 6)
add_edge!(g, 5, 7)
add_edge!(g, 3, 7)

append!(edge_list, [(3, 5), (1, 5), (5, 6), (5, 7), (3, 7)])
node_count += 3

arc_cost = 1e6 * ones(node_count, node_count)

travel_costs = zeros(edge_count + 5)
travel_costs[1:edge_count] = tn.free_flow_time
travel_costs[edge_count + 1] = 1
travel_costs[edge_count + 2] = 2
travel_costs[edge_count + 3] = 3
travel_costs[edge_count + 4] = 4
travel_costs[edge_count + 5] = 2
edge_count += 5

for (idx, (i,j)) in enumerate(edge_list)
    arc_cost[i, j] = travel_costs[idx]

    if j < first_non_zone_node
        arc_cost[i,j] = 1e6
    end
end

zone_outgoing_edges = [(i,j) for (i,j) in edge_list if i <=  od_pair_count]
zone_incoming_edges = [(i, j) for (i, j) in edge_list if j <= od_pair_count]
non_zone_edges = [(i, j) for (i, j) in edge_list if i > od_pair_count && j > od_pair_count]

removed_edges = [(3,4)]
# count_arcs_to_remove = Int(ceil(0.01 * length(edge_list)))
# removed_edges = StatsBase.sample(non_zone_edges, count_arcs_to_remove, replace=false)
removed_edge_idxs = [findall((x)->x==i, edge_list)[1] for i in removed_edges]
removed_edges_count = length(removed_edges)
y_val = zeros(removed_edges_count)

for edge in removed_edges
    rem_edge!(g, edge[1], edge[2])
end

in_neighbours = Dict()
out_neighbours = Dict()

for i in 1:node_count
    in_neighbours[i] = []
    out_neighbours[i] = []
end

for edge in edge_list
    push!(out_neighbours[edge[1]], edge[2])
    push!(in_neighbours[edge[2]], edge[1])
end


function get_edge_index(e)
    return findfirst(isequal(e), removed_edges)
end

nr_zone_outgoing_edges = [e for e in zone_outgoing_edges if !(e in removed_edges)]
nr_zone_incoming_edges = [e for e in zone_incoming_edges if !(e in removed_edges)]
nr_non_zone_edges = [e for e in non_zone_edges if !(e in removed_edges)]

r_zone_outgoing_edges = [e for e in zone_outgoing_edges if (e in removed_edges)]
r_zone_incoming_edges = [e for e in zone_incoming_edges if (e in removed_edges)]
r_non_zone_edges = [e for e in non_zone_edges if (e in removed_edges)]

model = Model(Gurobi.Optimizer) 

@variable(model, q[1:node_count, 1:od_pair_count])
@variable(model, s[1:removed_edges_count, 1:od_pair_count] >= 0)
@variable(model, r[1:od_pair_count, 1:od_pair_count])
@variable(model, t[1:od_pair_count] == 0)

@constraint(
    model, 
    potential_across_non_zone_arcs[(i,j) in nr_non_zone_edges, z in 1:od_pair_count],
    q[i,z] - q[j,z] <= arc_cost[i,j]
    )

@constraint(
    model, 
    potential_across_outgoing_arcs[(i,j) in nr_zone_outgoing_edges, z in 1:od_pair_count],
    r[i,z] - q[j,z] <= arc_cost[i,j]
    )

@constraint(
    model, 
    potential_across_incoming_arcs[(i,z) in nr_zone_incoming_edges],
    q[i,z] - t[z] <= arc_cost[i,z]
    )
       
@constraint(
    model, 
    potential_across_removed_non_zone_arcs[(i,j) in r_non_zone_edges, z in 1:od_pair_count],
    q[i,z] - q[j,z] - s[get_edge_index((i,j)),z] <= arc_cost[i,j]
    )

@constraint(
    model, 
    potential_across_removed_outgoing_arcs[(i,j) in r_zone_outgoing_edges, z in 1:od_pair_count],
    r[i,z] - q[j,z] - s[get_edge_index((i,j)),z] <= arc_cost[i,j]
    )

@constraint(
    model, 
    potential_across_removed_incoming_arcs[(i,z) in r_zone_incoming_edges],
    q[i,z] - t[z] - s[get_edge_index((i,z)),z] <= arc_cost[i,z]
    )

td_2_dest = [sum(demands[:, z]) for z in 1:od_pair_count] #Total travel demand to each destination

travel_cost_expression = @expression(
    model,
    sum(
        sum((r[i,z] - t[z]) * demands[i,z] for i in 1:od_pair_count) 
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

optimize!(model)


dists, parents = get_potentials_fw_mod(g, tn.free_flow_time, edge_list, tn.first_thru_node)

zdists = compute_zonal_distances(dists, parents, tn.free_flow_time, edge_list, od_pair_count)

potentials = compute_potentials(dists, parents, node_count, 
                                od_pair_count, tn.free_flow_time, edge_list, arc_cost, 
                                in_neighbours, out_neighbours)

s_vals = compute_s_vals(potentials, removed_edges, od_pair_count, arc_cost)
