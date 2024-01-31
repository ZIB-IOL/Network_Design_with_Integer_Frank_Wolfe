"""
This function updates the graphs by removing edges that are set to 0 and adding edges that are set to 1.
"""
function update_graph(o::Optimizer, idx, lb, ub)
    println("fixing index $(idx) to bounds $(lb), $(ub)")
    @assert lb == 0 || lb == 1 || isnothing(lb)
    @assert ub == 0 || ub == 1 || isnothing(ub)

    if isnothing(lb) && isnothing(ub)
        throw(
            ArgumentError(
                "Cannot set variable bounds to $(lb) and $(ub).",
            ),
        )
    end

    if !isnothing(lb) && !isnothing(ub)
        @assert lb <= ub
    end
    
    temp_rem_edges = []
    temo_add_edges = []

    if ub == 0.0
        tidx = idx - length(o.gparams.init_nodes) * (1 + o.gparams.od_pair_count)
        e = o.gparams.removed_edges[tidx]
        push!(temp_rem_edges, e)
        rem_edge!(o.graph, e[1], e[2])
        println("removing edge $(e)")
    end

    if lb == 1.0
        tidx = idx - length(o.gparams.init_nodes) * (1 + o.gparams.od_pair_count)
        e = o.gparams.removed_edges[tidx]
        push!(temo_add_edges, e)
        add_edge!(o.graph, e[1], e[2])
        println("adding edge $(e)")
    end

    return temp_rem_edges
end

function get_lmo_bounds(o)
    """
    This function obtains the upper and lower bounds from the SCIP.Optimizer model.
    """
    consLT_list =
        MOI.get(o, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.LessThan{Float64}}())
    consGT_list =
        MOI.get(o, MOI.ListOfConstraintIndices{MOI.VariableIndex,MOI.GreaterThan{Float64}}())

    upper_bounds = []
    lower_bounds = []
    for (i, cidx) in enumerate(consLT_list)
        if cidx.value in o.integer_variables
            push!(upper_bounds, MOI.get(o, MOI.ConstraintSet(), cidx).upper)
        end
    end
    for (i, cidx) in enumerate(consGT_list)
        if cidx.value in o.integer_variables
            push!(lower_bounds, MOI.get(o, MOI.ConstraintSet(), cidx).lower)
        end
    end

    return lower_bounds, upper_bounds
end

function modify_graph_with_bounds(o::Optimizer)
    """ 
    This function obtains the upper and lower bounds from the SCIP.Optimizer model. 
    It then updates the graph by removing edges that are set to 0 and adding edges that are set to 1.
    If the lower bound is different than the upper bound, the variable is added to the relaxed variables. 
    """
    
    lower_bounds, upper_bounds = get_lmo_bounds(o)

    @assert length(lower_bounds) == length(upper_bounds)

    fixed_variables = []
    relaxed_variables = []
    removed_edges = []
    for (i, lb, ub) in zip(o.integer_variables, lower_bounds, upper_bounds)
        if lb == ub
            push!(fixed_variables, i)
            eg = update_graph(o, i, lb, lb)
            push!(removed_edges, eg[1])
            #println("Fixed variable $i to $lb")
        elseif lb < ub
            push!(relaxed_variables, i)
        else
            throw(
                ArgumentError(
                    "Infeasible problem: lower bound is greater than upper bound for variable $i",
                ),
            )
        end
    end

    println("removed edges: $(removed_edges))")
    return fixed_variables, relaxed_variables, removed_edges
end

function update_binary_soln(o::Optimizer, v::Vector{Float64}, fixed_variables::Vector{}, relaxed_variables::Vector{})
    """
    This function updates the binary solution on the basis of the fixed and relaxed variables.
    It also verified the correctness of the solution.
    """
    lower_bounds, upper_bounds = get_lmo_bounds(o)

    for (i,iv) in enumerate(o.integer_variables)
        if iv in fixed_variables
            @assert o.curr_binary_soln[i] in [0,1]
            @assert o.curr_binary_soln[i] == upper_bounds[i] == lower_bounds[i]
        elseif iv in relaxed_variables
            removed_edge = o.gparams.removed_edges[i]
            removed_edge_idx = o.gparams.link_dic[removed_edge[1], removed_edge[2]]
            dst_arcs = o.gparams.od_pair_count * length(o.gparams.init_nodes)

            o.curr_binary_soln[i] = v[dst_arcs + removed_edge_idx] / o.gparams.MAX_FLOW
            #o.curr_solution[iv] = v[o.gparams.od_pair_count * length(o.gparams.init_nodes) + i] / o.gparams.MAX_FLOW
            #println("variable $(iv) is relaxed and has value $(o.curr_binary_soln[i])")
        end
    end
end


"""
The idea behind this function is to take a graph and optimizer. From the optimizer, read the upper and
lower bounds. If the bounds fix some of the edges then we either remove the edges or let them be. 
We also modify the input cost function to remove the elements corresponding to the removed edges. 
Then we solve the network flow problem over the reduced graph using this modified cost function. 
We map the returned solution onto a vector of length equal to the number of edges in the original graph. 
Finally, we add back the edges that have been removed and return the original solution. 
"""
function optimize_over_network_old(o::Optimizer, direction)
    edge_count = Graphs.ne(o.graph)
    println("_______________________________________")

    if length(direction) == 0
        var_count = edge_count + o.gparams.od_pair_count * edge_count
        println("vc:$(var_count)")     
        direction = collect(1:var_count)
    end

    println("initial edges in graph: $(Graphs.ne(o.graph))")
    ie = Graphs.ne(o.graph)
    println("ld: $(length(direction))")
    
    od_pair_count = o.gparams.od_pair_count
    
    function f(x)
        dest_edge_vars = x[1:(edge_count * od_pair_count)]
        edge_vars = x[(edge_count * od_pair_count + 1):(edge_count * od_pair_count + edge_count)]
        new_edge_vars = x[(edge_count * od_pair_count + edge_count + 1):end]
        return sum(edge_vars) + sum(o.gparams.removed_edges_costs .* new_edge_vars)
    end

    println("edges in graph: $(Graphs.ne(o.graph))")
    fixed_vars, relaxed_vars, removed_edges = modify_graph_with_bounds(o)
    println("edges in graph: $(Graphs.ne(o.graph))")
    removed_edge_idxs = [o.gparams.link_dic[eg[1],eg[2]] for eg in removed_edges]

    new_direction_end = copy(direction)[(ie * od_pair_count + 1):(ie * od_pair_count + ie)]
    new_direction_begin = copy(direction)[1:(ie * od_pair_count)]

    if length(removed_edges) > 0
        removed_edge_idxs = Int.(sort(removed_edge_idxs))
        println(removed_edge_idxs)
        deleteat!(new_direction_end, removed_edge_idxs)
        all_idx_to_remove = []

        for j in 1:(length(new_direction_begin) / ie)
            all_idx_to_remove = [all_idx_to_remove; (ie * (j - 1)) .+ removed_edge_idxs]
        end

        println(all_idx_to_remove)
        println(length(new_direction_begin))

        all_idx_to_remove = Int.(sort(all_idx_to_remove))
        deleteat!(new_direction_begin, all_idx_to_remove)
    end

    @assert length(new_direction_begin) == Graphs.ne(o.graph) * od_pair_count
    @assert length(new_direction_end) == Graphs.ne(o.graph)

    #println("lmofixed_variables:$(fixed_vars)")
    #println("lmorelaxed_variables:$(relaxed_vars)")
    #println("edges in graph $(length(collect(edges(o.graph))))")

    # for (i,iv) in enumerate(relaxed_vars)
    #     relx_edge = o.gparams.removed_edges[iv - var_count]
    #     println("edge modified $(relx_edge)")
    #     edge_idx = findall(x -> x == relx_edge, edge_list)[1]
    #     direction[o.gparams.od_pair_count * length(o.gparams.init_nodes) + edge_idx] += (o.gparams.removed_edges_costs[iv - var_count] / lmo.o.gparams.MAX_FLOW)
    #     #direction[iv - length(params.init_nodes)] += (o.gparams.removed_edges_costs[i] / o.gparams.MAX_FLOW)
    # end

    println("direction: $(round.(direction[end - 85:end], digits = 2))")
    dest_var_cnt = edge_count * o.gparams.od_pair_count
    #costs = direction[1:dest_var_cnt + edge_count]
    costs = [new_direction_begin; new_direction_end]
    v = all_or_nothing(costs, o.gparams, o.graph)
    
    orig_v = map_to_original(o, v)

    @assert all([orig_v[edge_count * o.gparams.od_pair_count + re] == 0.0 for re in removed_edge_idxs])
    @assert length(v) == Graphs.ne(o.graph) * o.gparams.od_pair_count + Graphs.ne(o.graph)

    zero_edges_zero = []
    for re in removed_edges
        edge_idx = o.gparams.link_dic[re[1], re[2]]
        edge_val = orig_v[ie * o.gparams.od_pair_count + edge_idx]
        push!(zero_edges_zero, edge_val == 0.0)
    end
    @assert all(zero_edges_zero)

    if v == MOI.INFEASIBLE
        o.curr_status = "Infeasible"
    else
        o.curr_status = "Solved"
    end

    term_st = MOI.get(o, MOI.TerminationStatus())

    if term_st ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.SLOW_PROGRESS)
        @error "Unexpected termination: $term_st"
        return o.curr_binary_soln
    end

    for eg in removed_edges
        #println("Adding edge $(eg[1]) -> $(eg[2])")
        t = Graphs.add_edge!(o.graph, eg[1], eg[2])
        #println("edge added: $(t)")
    end

    #update the binary solution for relaxed variables
    update_binary_soln(o, orig_v, fixed_vars, relaxed_vars)
    o.curr_solution = [orig_v; o.curr_binary_soln]

    println("v eval $(sum(v .* costs))")
    println("cep f(v): $(f([orig_v; o.curr_binary_soln]))")


    println("final edges in graph $(Graphs.ne(o.graph))")
    @assert Graphs.ne(o.graph) == ie
    return [orig_v; o.curr_binary_soln] 
end

function map_to_original(o::Optimizer, v::Vector{Float64})
    """
    Maps the solution of the all-or-nothing problem to a solution of the original problem
    """
    @assert v != MOI.INFEASIBLE

    orig_edge_list = collect(zip(o.gparams.init_nodes, o.gparams.term_nodes))
    curr_edge_list = [(eg.src, eg.dst) for eg in edges(o.graph)]

    orig_v = zeros(length(orig_edge_list) + o.gparams.od_pair_count * length(orig_edge_list))

    orig_arc_start = o.gparams.od_pair_count * length(orig_edge_list)
    curr_arc_start = o.gparams.od_pair_count * length(curr_edge_list)

    for (i,eg) in enumerate(curr_edge_list)
        edge_idx = findfirst(x -> x == eg, orig_edge_list)[1]
        #println("searching for edge $eg")
        #println("edi $(edge_idx)")
        orig_v[orig_arc_start + edge_idx] = v[curr_arc_start + i]
    end

    for dst in 1:o.gparams.od_pair_count
        for (i,eg) in enumerate(curr_edge_list)
            edge_idx = findfirst(x -> x == eg, orig_edge_list)[1]
            #println("searching for edge $eg")
            #println("edi $(edge_idx)")
            orig_v[(dst - 1) * length(orig_edge_list) + edge_idx] = v[(dst - 1) * length(curr_edge_list) + i]
        end
    end

    @assert all(orig_v .>= 0)

    for edge in orig_edge_list
        agg_x_flow = orig_v[length(orig_edge_list) * o.gparams.od_pair_count + o.gparams.link_dic[edge[1],edge[2]]]
        sum_x_flow = sum(orig_v[o.gparams.link_dic[edge[1],edge[2]] + (j-1) * length(orig_edge_list)] for j in 1:o.gparams.od_pair_count)
        
        @assert agg_x_flow ≈ sum_x_flow
    end

    return orig_v
end