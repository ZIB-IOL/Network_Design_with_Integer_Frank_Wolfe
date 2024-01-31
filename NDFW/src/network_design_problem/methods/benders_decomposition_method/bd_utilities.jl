"""
Get value of the dual variable s depending upon y
When y = 1, s is set to 0
when y = 0, s is set to the cost of the removed edge
"""
function get_dual_var_s(o, y_vals, direc)
    cec = length(o.gparams.init_nodes)
    rec = length(o.gparams.removed_edges)
    opc = o.gparams.od_pair_count
    edge_costs = direc
    dual_var_s = Dict()

    for i in eachindex(y_vals)
        yv = y_vals[i]
        cor_edge = o.gparams.removed_edges[i]
        cor_edge_idx = o.gparams.link_dic[cor_edge[1], cor_edge[2]]
        if yv == 1.0
            dual_var_s[cor_edge_idx] = zeros(opc)
        elseif yv == 0.0
            #@info edge_costs[cor_edge_idx]
            dual_var_s[cor_edge_idx] = edge_costs[cor_edge_idx] * ones(opc)

            if cor_edge[2] <= o.gparams.first_thru_node
                dual_var_s[cor_edge_idx] = 1e4 * ones(opc)
            end


            #dual_var_s[cor_edge_idx] = ones(opc)
            # for z in 1:opc
            #     dual_var_s[cor_edge_idx] = max(potentials[(cor_edge[1], z)] - potentials[(cor_edge[2], z)] - edge_costs[cor_edge_idx], 0.0)
            # end
        else
            println("y_vals[$i] = $yv")
            throw(ArgumentError("Incorrect value for y."))
        end
    end

    return dual_var_s
end


# function get_vars_by_state(o::Optimizer)
#     """ 
#     This function returns the integer variables in two categories: fixed and relaxed.
#     """

#     lower_bounds, upper_bounds = get_lmo_bounds(o)

#     @assert length(lower_bounds) == length(upper_bounds)

#     fixed_variables = []
#     relaxed_variables = []

#     for (i, lb, ub) in zip(o.integer_variables, lower_bounds, upper_bounds)
#         if lb == ub
#             push!(fixed_variables, i)
#             #println("Fixed variable $i to $lb")
#         elseif lb < ub
#             push!(relaxed_variables, i)
#         else
#             throw(
#                 ArgumentError(
#                     "Infeasible problem: lower bound is greater than upper bound for variable $i",
#                 ),
#             )
#         end
#     end

#     return fixed_variables, relaxed_variables
# end

"""
This function takes in a solution to the sub problem and returns the upper bound on the objective value.
INPUT
o: Boscia object
v: solution to the subproblem
costs: costs of the edges
bin_var_costs: costs of the binary variables
y_vals: values of the binary variables
"""
function get_obj_upper_bound(gparams, orig_r, orig_t, orig_s, costs, bin_var_costs, y_vals)
    #edge_count = Graphs.ne(o.graph)
    demands = gparams.travel_demand
    total_scenarios = gparams.total_scenarios
    od_pair_count = gparams.od_pair_count

    #travel_costs = sum([costs[i] * v[i] for i in edge_count* od_pair_count + 1:edge_count * (od_pair_count + 1)]) 

    demand_costs = 0.0
    for scen in 1:total_scenarios
        demand_costs += sum(
            sum((orig_r[scen][i, z] - orig_t[scen][z]) * demands[scen][i, z] for i in 1:od_pair_count)
            for z in 1:od_pair_count
        )
    end

    td_2_dest = [maximum([sum(gparams.travel_demand[s][:, z]) for s in 1:gparams.total_scenarios]) for z in 1:gparams.od_pair_count]

    dual_costs = 0.0
    for scen in 1:total_scenarios
        terms_vals = sum([td_2_dest[z] * orig_s[scen][e, z] * y_vals[e] for e in eachindex(y_vals) for z in 1:gparams.od_pair_count])
        dual_costs += terms_vals
    end

    arc_costs = sum([bin_var_costs[i] * y_vals[i] for i in eachindex(bin_var_costs)])
    upperbound = demand_costs + arc_costs - dual_costs

    return upperbound
end


"""
Get the edges that cross over from the component connected to the 
origin to the component NOT connected to the origin.  
"""
function get_bridge_edges(gparams, state, src, costs, orig_edge_list)
    distances = copy(state.dists)

    #We fix the costs of edges going into the zones because by default it is Inf
    for zone in 1:(gparams.first_thru_node-1)
        if state.parents[zone] == 0
            continue
        end
        edge = (state.parents[zone], zone)
        edge_idx = findall(x -> x == edge, orig_edge_list)[1]

        distances[zone] = distances[state.parents[zone]] + costs[edge_idx]
    end

    #nodes whose distance to the origin is finite
    o_graph = collect(1:length(distances))[distances.<Inf]
    #nodes whose distance to the origin is infinite
    d_graph = collect(1:length(distances))[distances.==Inf]

    #edges that cross over from the origin component (finite distance) 
    #to the destination component (infinite distance)
    bridge_edges = []

    for e in gparams.removed_edges
        if e[1] in o_graph && e[2] in d_graph
            push!(bridge_edges, e)
        end
    end

    return bridge_edges
end

"""
Get the edges that cross over from the component connected to the 
origin to the component NOT connected to the origin.  
Get these edges corresponding ot multiple steps
"""
function get_multi_step_bridge_edges(gparams, graph, state, src, costs, orig_edge_list)
    distances = copy(state.dists)
    bridge_edges = Dict()

    #We fix the costs of edges going into the zones
    for zone in 1:(gparams.first_thru_node-1)
        if state.parents[zone] == 0
            continue
        end
        edge = (state.parents[zone], zone)
        edge_idx = findall(x -> x == edge, orig_edge_list)[1]
        distances[zone] = distances[state.parents[zone]] + costs[edge_idx]
    end

    o_graph = collect(1:length(distances))[distances.<Inf]
    d_graph = collect(1:length(distances))[distances.==Inf]

    mod_graph = copy(graph)

    nodes_added = Dict()
    nodes_added_all = []
    for step in 1:length(orig_edge_list)
        bridge_edges[step] = []
        nodes_added[step] = []

        for e in gparams.removed_edges
            if step == 1
                if e[1] in o_graph && e[2] in d_graph
                    push!(bridge_edges[step], e)
                    push!(nodes_added[step], e[2])
                    push!(nodes_added_all, e[2])
                    add_edge!(mod_graph, e[1], e[2])
                end
            else
                comps = Graphs.strongly_connected_components(mod_graph)
                main_comp = [i for i in comps if src in i][1]

                if e[1] in main_comp && e[2] in !(e[2] in main_comp)
                    push!(bridge_edges[step], e)
                    push!(nodes_added[step], e[2])

                    if e[2] == dst
                        break
                    end
                end
            end
        end

        if length(bridge_edges[step]) == 0
            delete!(bridge_edges, step)
            delete!(nodes_added, step)
            break #TODO Verify that this is enough. 
        end
    end

    @info bridge_edges
    @info nodes_added

    return bridge_edges
end


# """
# Restore the graph to its original state by adding back the removed edges
# """
# function restore_graph(o::Optimizer, removed_edges)
#     for eg in removed_edges
#         #println("Adding edge $(eg[1]) -> $(eg[2])")
#         t = Graphs.add_edge!(o.graph, eg[1], eg[2])
#         #println("edge added: $(t)")
#     end
# end

# """
# Obtain new directions after removing the costs of the removed edges
# """
# function obtain_new_directions(o::Optimizer, direction, actually_removed_edges, ie)
#     removed_edges = actually_removed_edges
#     od_pair_count = o.gparams.od_pair_count
#     removed_edge_idxs = [o.gparams.link_dic[eg[1], eg[2]] for eg in removed_edges]

#     new_direction_end = copy(direction)[(ie*od_pair_count+1):(ie*od_pair_count+ie)]
#     new_direction_begin = copy(direction)[1:(ie*od_pair_count)]

#     if length(removed_edges) > 0
#         removed_edge_idxs = Int.(sort(removed_edge_idxs))
#         deleteat!(new_direction_end, removed_edge_idxs)
#         all_idx_to_remove = []

#         for j in 1:(length(new_direction_begin)/ie)
#             all_idx_to_remove = [all_idx_to_remove; (ie * (j - 1)) .+ removed_edge_idxs]
#         end


#         all_idx_to_remove = Int.(sort(all_idx_to_remove))
#         deleteat!(new_direction_begin, all_idx_to_remove)
#     end

#     return new_direction_begin, new_direction_end
# end

# """
# This function use the dual variables for the source constraint and the distance from the source
# to each node to calculate the potentials. 
# """
# function get_potentials(o, r, rev_pot, parents)
#     no_vertices = Graphs.nv(o.graph)
#     q = Dict()
#     ver_q = Dict()
#     for src in 1:(o.gparams.first_thru_node-1)
#         for dst in 1:(o.gparams.first_thru_node-1)
#             i = dst
#             dist_so_far = r[src][dst]
#             q[(i, dst)] = dist_so_far - rev_pot[src][i]

#             if haskey(ver_q, (i, dst))
#                 push!(ver_q[(i, dst)], q[(i, dst)])
#             else
#                 ver_q[(i, dst)] = [q[(i, dst)]]
#             end

#             while i != src
#                 i = parents[src][i]
#                 q[(i, dst)] = dist_so_far - rev_pot[src][i]

#                 if haskey(ver_q, (i, dst))
#                     push!(ver_q[(i, dst)], q[(i, dst)])
#                 else
#                     ver_q[(i, dst)] = [q[(i, dst)]]
#                 end

#             end
#         end
#     end

#     return q, ver_q
# end

# """
# This function use the dual variables for the source constraint and the distance from the source
# to each node to calculate the potentials. 
# """
# function get_potentials_fw(o, travel_time_all)
#     edge_count = Graphs.ne(o.graph)
#     od_pair_count = o.gparams.od_pair_count

#     edge_list = [(edge.src, edge.dst) for edge in collect(edges(o.graph))]
#     init_nodes = [i for (i, j) in edge_list]
#     term_nodes = [j for (i, j) in edge_list]
#     new_link_dic = sparse(init_nodes, term_nodes, collect(1:edge_count))


#     travel_time = travel_time_all[edge_count*od_pair_count+1:edge_count*(od_pair_count+1)]     #take only the aggregated components

#     state = TA_floyd_warshall_shortest_paths(o.graph, travel_time, new_link_dic, o.gparams.first_thru_node)
#     return state.dists
# end

"""
Normalize the cut by reducing the right hand side coefficient
"""
function normalize_cut(c, normalizing_consts)
    original_coeff = c[1][2].lower
    original_coeff_list[c] = original_coeff

    if length(normalizing_consts) == 1
        normalizing_const = minimum(normalizing_consts)
        MOI.set(bopt, MOI.ConstraintSet(), c, MOI.GreaterThan(original_coeff - normalizing_const))
    else
        if minimum(normalizing_consts) == normalizing_const
            original_coeff = original_coeff_list[c]
            MOI.set(bopt, MOI.ConstraintSet(), c, MOI.GreaterThan(original_coeff - normalizing_const))
        elseif minimum(normalizing_consts) < normalizing_const
            normalizing_const = minimum(normalizing_consts)
            for (idx, oc) in enumerate(optimality_cuts) 
                original_coeff = original_coeff_list[oc]
                MOI.set(bopt, MOI.ConstraintSet(), oc, MOI.GreaterThan(original_coeff - normalizing_const))
            end
            original_coeff =original_coeff_list[c]
            MOI.set(bopt, MOI.ConstraintSet(), c, MOI.GreaterThan(original_coeff - normalizing_const))
        end
    end
end

function check_for_optimality(lowerbound, gparams, orig_r, orig_t, orig_s, costs, bin_var_costs, y_vals, eps)
    #@info "Current solution: $(r), $(t), $(s)"
    upperbound = get_obj_upper_bound(gparams, orig_r, orig_t, orig_s, costs, bin_var_costs, y_vals)

    #potentials = get_potentials_fw(o, costs)
    #s = get_dual_var_s(o, y_vals, edge_costs)

    @info "gap $(lowerbound), $(upperbound), $(upperbound - lowerbound), $((upperbound - lowerbound)/abs(lowerbound))"

    if (upperbound - lowerbound)/abs(lowerbound) < eps         
        return true
    end

    return false
end