"""
Network Design LMO
"""
mutable struct BNDesignLMO <: Boscia.BoundedLinearMinimizationOracle
    n::Int64
    graph
    int_vars::Vector{Int}
    bounds::Boscia.IntegerBounds
    gparams::GraphParams
    use_bigm::Bool
    use_adaptive::Bool
    use_reverse::Bool
    solving_time::Float64
    distmx::Matrix{Float64}
end

BNDesignLMO(n, graph, int_vars, bounds, gparams, use_bigm, use_adaptive, use_reverse) = BNDesignLMO(n, graph, int_vars, bounds, gparams, use_bigm, use_adaptive, use_reverse, 0.0, zeros(Float64, Graphs.nv(graph), Graphs.nv(graph)))

function Boscia.compute_extreme_point(ndlmo::BNDesignLMO, direction; kwargs...)
    graph = ndlmo.graph
    gparams = ndlmo.gparams
    int_vars = ndlmo.int_vars

    ie = Graphs.ne(graph)
    od_pair_count = gparams.od_pair_count
    total_scenarios = gparams.total_scenarios
    use_bigm = ndlmo.use_bigm
    use_adaptive = ndlmo.use_adaptive
    use_reverse = ndlmo.use_reverse
    distmx = ndlmo.distmx

    total_scenario_vars = od_pair_count * ie + ie

    bin_var_costs = direction[int_vars]
    costs = direction[1:total_scenarios * total_scenario_vars]

    lower_bounds = zeros(length(int_vars))
    upper_bounds = ones(length(int_vars))

    for i in eachindex(int_vars)
        lower_bounds[i] = ndlmo.bounds[int_vars[i], :greaterthan]
        upper_bounds[i] = ndlmo.bounds[int_vars[i], :lessthan]
    end

    # for v in ndlmo.int_vars
    #     println("$(ndlmo.bounds[v,:greaterthan]) <= x[$(v)] <= $(ndlmo.bounds[v,:lessthan])")
    # end


    # Take the vector of costs and set the positive ones to 0 and negative ones to 1.
    # This is the minimizer. 
    #temp_results_y = 1.0 * (min.(0, bin_var_costs) .!= 0)
    temp_results_y = zeros(length(bin_var_costs))
    adaptive_y = ones(length(bin_var_costs))

    #We update the computed solution to account for the fixed variables. 
    #Since there are not constraints in the problem all this means is that 
    #we set the fixed variables to their fixed values.
    for i in eachindex(gparams.removed_edges)
        if lower_bounds[i] == upper_bounds[i]
            if upper_bounds[i] == 0.0
                temp_results_y[i] = 0.0
                adaptive_y[i] = 0.0
            elseif lower_bounds[i] == 1.0
                temp_results_y[i] = 1.0
                adaptive_y[i] = 1.0
            else
                throw(ArgumentError("Incorrect bounds for fixed variable."))
            end
        elseif lower_bounds[i] < upper_bounds[i]
            if bin_var_costs[i] < 0
                temp_results_y[i] = 1.0
            else
                temp_results_y[i] = 0.0
            end
        else
            throw(ArgumentError("Variable not in fixed or relaxed set."))
        end
    end


    if use_adaptive == true
        orig_v = get_adaptive_flow_extreme_point(adaptive_y, costs, graph, ie, gparams, use_bigm, use_reverse, distmx)
    elseif use_adaptive == false
        orig_v = get_static_flow_extreme_point(costs, graph, ie, gparams, use_bigm, use_reverse, distmx)
    else
        throw(ArgumentError("Invalid use_adaptive argument $(use_adaptive)"))
    end


    #This has to be after mapping v to the original vector
    #restore_graph(o, actually_removed_edges)


    #@assert Boscia.is_linear_feasible(o.sopt, [orig_v; temp_results_y])

    # if v == MOI.INFEASIBLE
    #     o.curr_status = "Infeasible"
    # else
    #     o.curr_status = "Solved"
    # end

    #term_st = MOI.get(o, MOI.TerminationStatus())

    # if term_st ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.SLOW_PROGRESS)
    #     @error "Unexpected termination: $term_st"
    #     return o.curr_binary_soln
    # end

    curr_solution = [orig_v; temp_results_y]

    @assert Graphs.ne(graph) == ie

    return curr_solution
end

function modify_graph(graph, y_vals, removed_edges)
    """ 
    This function obtains the upper and lower bounds from the SCIP.Optimizer model. 
    It then updates the graph by removing edges that are set to 0 and adding edges that are set to 1.
    If the lower bound is different than the upper bound, the variable is added to the relaxed variables. 
    """

    for i in eachindex(y_vals)
        @assert y_vals[i] == 0.0 || y_vals[i] == 1.0
    end

    currently_removed_edges = []
    for (i, yv) in enumerate(y_vals)
        if yv == 0.0
            edge = removed_edges[i]
            rem_edge!(graph, edge[1], edge[2])
            push!(currently_removed_edges, edge)
        elseif yv == 1.0
            continue
        else
            throw(
                ArgumentError(
                    "incorrect value of y",
                ),
            )
        end
    end

    return graph, currently_removed_edges
end

function modify_cost_vector(direction, gparams, curr_removed_edges, ie)

    if length(curr_removed_edges) == 0
        return direction
    end

    od_pair_count = gparams.od_pair_count
    removed_edge_idxs = [gparams.link_dic[eg[1], eg[2]] for eg in curr_removed_edges]
    total_scenarios = gparams.total_scenarios

    new_direction = copy(direction)

    ind_idxs_to_remove = []
    for scen in 1:total_scenarios
        first_var = (scen - 1) * ie * (od_pair_count + 1)
        for j in 1:od_pair_count
            ind_idxs_to_remove = [ind_idxs_to_remove; first_var .+ (ie * (j - 1)) .+ removed_edge_idxs]
        end
    end

    agg_removed_edge_idxs = []
    for scen in 1:total_scenarios
        first_var = (scen - 1) * ie * (od_pair_count + 1)
        agg_removed_edge_idxs = [agg_removed_edge_idxs; first_var .+ ie * od_pair_count .+ removed_edge_idxs]
    end


    removed_edge_idxs = [ind_idxs_to_remove; agg_removed_edge_idxs]

    if length(removed_edge_idxs) > 0
        removed_edge_idxs = Int.(sort(removed_edge_idxs))
        deleteat!(new_direction, removed_edge_idxs)
    end

    return new_direction
end

function map_to_original(v::Vector{Float64}, graph, gparams, total_scenarios)
    """
    Maps the solution of the all-or-nothing problem to a solution of the original problem
    """
    @assert v != MOI.INFEASIBLE

    od_pair_count = gparams.od_pair_count

    orig_edge_list = collect(zip(gparams.init_nodes, gparams.term_nodes))
    curr_edge_list = [(eg.src, eg.dst) for eg in edges(graph)]

    total_scenario_variables = length(orig_edge_list) + gparams.od_pair_count * length(orig_edge_list)
    orig_v = zeros(total_scenarios * total_scenario_variables)

    for scen in 1:total_scenarios

        first_var = (scen - 1) * total_scenario_variables
        curr_first_var = (scen - 1) * (length(curr_edge_list) * (od_pair_count + 1))

        
        orig_arc_start = od_pair_count * length(orig_edge_list)
        curr_arc_start = od_pair_count * length(curr_edge_list)

        for (i, eg) in enumerate(curr_edge_list)
            edge_idx = findfirst(x -> x == eg, orig_edge_list)[1]
            orig_v[first_var + orig_arc_start + edge_idx] = v[curr_first_var + curr_arc_start + i]
        end

        for dst in 1:gparams.od_pair_count
            for (i, eg) in enumerate(curr_edge_list)
                edge_idx = findfirst(x -> x == eg, orig_edge_list)[1]
                orig_v[first_var+(dst-1)*length(orig_edge_list)+edge_idx] = v[curr_first_var + (dst - 1)*length(curr_edge_list)+i]
            end
        end

    end
    @assert all(orig_v .>= 0)

    oel_count = length(orig_edge_list)
    for edge in orig_edge_list
        agg_x_flow = orig_v[oel_count*gparams.od_pair_count+gparams.link_dic[edge[1], edge[2]]]
        sum_x_flow = sum(orig_v[gparams.link_dic[edge[1], edge[2]]+(j-1)*oel_count] for j in 1:gparams.od_pair_count)
        @assert agg_x_flow ≈ sum_x_flow
    end

    return orig_v
end

function get_static_flow_extreme_point(direction, graph, ie, gparams, use_bigm, use_reverse, distmx)
    od_pair_count = gparams.od_pair_count
    total_scenarios = gparams.total_scenarios
    total_scenario_vars = ie * (od_pair_count + 1)
    orig_edge_list = [(gparams.init_nodes[i], gparams.term_nodes[i]) for i in 1:length(gparams.init_nodes)]

    orig_v = zeros(total_scenario_vars * total_scenarios)

    temp_graph = copy(graph)

    reverse_temp_graph = create_reverse_graph(temp_graph)

    for scen = 1:total_scenarios
        first_var = (scen - 1) * total_scenario_vars
        costs = direction[(first_var+1):(first_var+(od_pair_count*ie)+ie)]

        #v, _, _ = all_or_nothing_stochastic_detailed(costs, gparams, graph, graph, scen)
        if use_bigm == false
            if use_reverse == false
                v = all_or_nothing_stochastic_separated(costs, gparams, temp_graph, orig_edge_list, scen, distmx)
            elseif use_reverse == true
                v = all_or_nothing_stochastic_reversed(costs, gparams, reverse_temp_graph, orig_edge_list, scen, distmx)
            else
                throw(ArgumentError("Invalid use_reverse argument $(use_reverse)"))
            end
        else
            v = all_or_nothing_stochastic_common(costs, gparams, temp_graph, temp_graph, scen, distmx)
        end

        # if v == MOI.INFEASIBLE
        #     curr_status = "Infeasible"
        #     break
        # else
        #     curr_status = "Solved"
        # end
        orig_v[(first_var+1):(first_var+(od_pair_count*ie)+ie)] .= v
    end


    return orig_v
end

function get_adaptive_flow_extreme_point(temp_results_y, direction, graph, ie, gparams, use_bigm, use_reverse, distmx)
    od_pair_count = gparams.od_pair_count
    total_scenarios = gparams.total_scenarios
    total_scenario_vars = ie * (od_pair_count + 1)
    orig_edge_list = [(gparams.init_nodes[i], gparams.term_nodes[i]) for i in 1:length(gparams.init_nodes)]

    orig_v = zeros(total_scenario_vars * total_scenarios)

    temp_graph = copy(graph)

    modified_graph, curr_removed_edges = modify_graph(temp_graph, temp_results_y, gparams.removed_edges)

    curr_removed_count = size(curr_removed_edges)[1]
    current_edges = ie - curr_removed_count
    modified_total_scenario_vars = od_pair_count * current_edges + current_edges
    modified_costs = modify_cost_vector(direction, gparams, curr_removed_edges, ie)

    curr_edge_list = copy(orig_edge_list)
    removed_edge_idxs = [gparams.link_dic[eg[1], eg[2]] for eg in curr_removed_edges]

    if length(removed_edge_idxs) > 0
        removed_edge_idxs = Int.(sort(removed_edge_idxs))
        deleteat!(curr_edge_list, removed_edge_idxs)
    end

    reverse_temp_graph = create_reverse_graph(temp_graph)
    
    @assert length(modified_costs) == total_scenarios * modified_total_scenario_vars
    
    orig_v = zeros(modified_total_scenario_vars * total_scenarios)

    @assert Graphs.ne(temp_graph) == length(curr_edge_list)

    #new_direction_begin, new_direction_end = obtain_new_directions(o, res_direction, actually_removed_edges, ie)
    for scen = 1:total_scenarios
        first_var = (scen - 1) * modified_total_scenario_vars
        costs = modified_costs[(first_var+1):(first_var+(od_pair_count*current_edges)+current_edges)]

        #v, _, _ = all_or_nothing_stochastic_detailed(costs, gparams, graph, graph, scen)
        if use_bigm == false
            if use_reverse == false
                v = all_or_nothing_stochastic_separated(costs, gparams, temp_graph, curr_edge_list, scen, distmx)
            elseif use_reverse == true
                v = all_or_nothing_stochastic_reversed(costs, gparams, reverse_temp_graph, curr_edge_list, scen, distmx)
            else
                throw(ArgumentError("Invalid use_reverse argument $(use_reverse)"))
            end
        elseif use_bigm == true
            v = all_or_nothing_stochastic_common(costs, gparams, temp_graph, graph, scen, distmx)
        else
            throw(ArgumentError("Invaid use_bigm value $(use_bigm)"))
        end

        # if v == MOI.INFEASIBLE
        #     curr_status = "Infeasible"
        #     break
        # else
        #     curr_status = "Solved"
        # end
        orig_v[(first_var+1):(first_var+(od_pair_count*current_edges)+current_edges)] .= v
    end


    final_v = map_to_original(orig_v, temp_graph, gparams, total_scenarios)
    #
    return final_v
end

function check_feasible_flow(temp_results_y, direction, graph, ie, gparams, use_bigm, use_reverse, distmx)
    od_pair_count = gparams.od_pair_count
    total_scenarios = gparams.total_scenarios
    total_scenario_vars = ie * (od_pair_count + 1)
    orig_edge_list = [(gparams.init_nodes[i], gparams.term_nodes[i]) for i in 1:length(gparams.init_nodes)]

    orig_v = zeros(total_scenario_vars * total_scenarios)

    temp_graph = copy(graph)

    modified_graph, curr_removed_edges = modify_graph(temp_graph, temp_results_y, gparams.removed_edges)
    
    curr_edge_list = copy(orig_edge_list)
    removed_edge_idxs = [gparams.link_dic[eg[1], eg[2]] for eg in curr_removed_edges]

    if length(removed_edge_idxs) > 0
        removed_edge_idxs = Int.(sort(removed_edge_idxs))
        deleteat!(curr_edge_list, removed_edge_idxs)
    end

    
    curr_removed_count = size(curr_removed_edges)[1]
    current_edges = ie - curr_removed_count
    modified_total_scenario_vars = od_pair_count * current_edges + current_edges
    modified_costs = modify_cost_vector(direction, gparams, curr_removed_edges, ie)

    reverse_temp_graph = create_reverse_graph(temp_graph)


    @assert length(modified_costs) == total_scenarios * modified_total_scenario_vars

    @assert Graphs.ne(temp_graph) == length(curr_edge_list)

    orig_v = zeros(modified_total_scenario_vars * total_scenarios)

    #new_direction_begin, new_direction_end = obtain_new_directions(o, res_direction, actually_removed_edges, ie)
    for scen = 1:total_scenarios
        first_var = (scen - 1) * modified_total_scenario_vars
        costs = modified_costs[(first_var+1):(first_var+(od_pair_count*current_edges)+current_edges)]

        #v, _, _ = all_or_nothing_stochastic_detailed(costs, gparams, graph, graph, scen)
        if use_bigm == false
            if use_reverse == false
                v = all_or_nothing_stochastic_separated(costs, gparams, temp_graph, curr_edge_list, scen, distmx)
            elseif use_reverse == true
                v = all_or_nothing_stochastic_reversed(costs, gparams, reverse_temp_graph, curr_edge_list, scen, distmx)
            else
                throw(ArgumentError("Invalid use_reverse argument $(use_reverse)"))
            end
        elseif use_bigm == true
            v = all_or_nothing_stochastic_common(costs, gparams, temp_graph, graph, scen, distmx)
        else
            throw(ArgumentError("Invaid use_bigm value $(use_bigm)"))
        end

        if v == MOI.INFEASIBLE
            return MOI.INFEASIBLE
        end

        orig_v[(first_var+1):(first_var+(od_pair_count*current_edges)+current_edges)] .= v
    end


    final_v = map_to_original(orig_v, temp_graph, gparams, total_scenarios)

    return final_v
end

function Boscia.build_global_bounds(ndlmo::BNDesignLMO, integer_variables)
    global_bounds = Boscia.IntegerBounds()
    for i in 1:ndlmo.n
        if i in integer_variables
            push!(global_bounds, (i, ndlmo.bounds[i, :lessthan]), :lessthan)
            push!(global_bounds, (i, ndlmo.bounds[i, :greaterthan]), :greaterthan)
        end
    end
    return global_bounds
end

## Read information from problem
function Boscia.get_list_of_variables(ndlmo::BNDesignLMO)
    return ndlmo.n, collect(1:ndlmo.n)
end

# Get list of integer variables, respectively.
function Boscia.get_integer_variables(ndlmo::BNDesignLMO)
    return ndlmo.int_vars
end

function Boscia.get_int_var(ndlmo::BNDesignLMO, cidx)
    return cidx
end

function Boscia.get_lower_bound_list(ndlmo::BNDesignLMO)
    return keys(ndlmo.bounds.lower_bounds)
end

function Boscia.get_upper_bound_list(ndlmo::BNDesignLMO)
    return keys(ndlmo.bounds.upper_bounds)
end

function Boscia.get_bound(ndlmo::BNDesignLMO, c_idx, sense::Symbol)
    @assert sense == :lessthan || sense == :greaterthan
    return ndlmo[c_idx, sense]
end

## Changing the bounds constraints.
function Boscia.set_bound!(ndlmo::BNDesignLMO, c_idx, value, sense::Symbol)
    if sense == :greaterthan
        ndlmo.bounds.lower_bounds[c_idx] = value
    elseif sense == :lessthan
        ndlmo.bounds.upper_bounds[c_idx] = value
    else
        error("Allowed values for sense are :lessthan and :greaterthan.")
    end
end

function Boscia.delete_bounds!(ndlmo::BNDesignLMO, cons_delete)
    # For the cube this shouldn't happen! Otherwise we get unbounded!
    if !isempty(cons_delete)
        error("Trying to delete bounds of the cube!")
    end
end

function Boscia.add_bound_constraint!(ndlmo::BNDesignLMO, key, value, sense::Symbol)
    # Should not be necessary
    return error("Trying to add bound constraints of the cube!")
end

## Checks

function Boscia.is_constraint_on_int_var(ndlmo::BNDesignLMO, c_idx, int_vars)
    return c_idx in int_vars
end

function Boscia.is_bound_in(ndlmo::BNDesignLMO, c_idx, bounds)
    return haskey(bounds, c_idx)
end

function Boscia.is_linear_feasible(ndlmo::BNDesignLMO, v::AbstractVector)
    for i in ndlmo.int_vars
        if !(
            ndlmo.bounds[i, :greaterthan] ≤ v[i] + 1e-6 ||
            !(v[i] - 1e-6 ≤ ndlmo.bounds[i, :lessthan])
        )
            @debug(
                "Vertex entry: $(v[i]) Lower bound: $(ndlmo.bounds[i, :greaterthan]) Upper bound: $(ndlmo.bounds[i, :lessthan]))"
            )
            return false
        end
    end
    return true
end

function Boscia.has_integer_constraint(ndlmo::BNDesignLMO, idx)
    return idx in ndlmo.int_vars
end


###################### Optional
## Safety Functions

function Boscia.build_LMO_correct(ndlmo::BNDesignLMO, node_bounds)
    for key in keys(node_bounds.lower_bounds)
        if !haskey(ndlmo.bounds, (key, :greaterthan)) ||
           ndlmo.bounds[key, :greaterthan] != node_bounds[key, :greaterthan]
            return false
        end
    end
    for key in keys(node_bounds.upper_bounds)
        if !haskey(ndlmo.bounds, (key, :lessthan)) ||
           ndlmo.bounds[key, :lessthan] != node_bounds[key, :lessthan]
            return false
        end
    end
    return true
end

function Boscia.check_feasibility(ndlmo::BNDesignLMO)
    #println("Feasiblity test started")
    for i in ndlmo.int_vars
        if !haskey(ndlmo.bounds, (i, :greaterthan)) || !haskey(ndlmo.bounds, (i, :lessthan))
            return Boscia.UNBOUNDED
        end
        if ndlmo.bounds[i, :greaterthan] > ndlmo.bounds[i, :lessthan]
            return Boscia.INFEASIBLE
        end
    end

    int_vars = ndlmo.int_vars
    gparams = ndlmo.gparams
    ie = length(gparams.init_nodes)
    use_bigm = ndlmo.use_bigm
    use_reverse = ndlmo.use_reverse
    od_pair_count = gparams.od_pair_count
    total_scenarios = gparams.total_scenarios

    lower_bounds = zeros(length(int_vars))
    upper_bounds = ones(length(int_vars))

    for i in eachindex(int_vars)
        lower_bounds[i] = ndlmo.bounds[int_vars[i], :greaterthan]
        upper_bounds[i] = ndlmo.bounds[int_vars[i], :lessthan]
    end

    temp_results_y = ones(length(int_vars))

    for i in eachindex(gparams.removed_edges)
        if lower_bounds[i] == upper_bounds[i]
            if upper_bounds[i] == 0.0
                temp_results_y[i] = 0.0
            elseif lower_bounds[i] == 1.0
                temp_results_y[i] = 1.0
            else
                throw(ArgumentError("Incorrect bounds for fixed variable."))
            end
        end
    end

    direction = ones(total_scenarios  * (od_pair_count * ie + ie))
    
    temp_graph = copy(ndlmo.graph)
    if ndlmo.use_adaptive == true
        v = check_feasible_flow(temp_results_y, direction, temp_graph, ie, gparams, use_bigm, use_reverse, ndlmo.distmx)
    elseif ndlmo.use_adaptive == false
        all_edges_added = ones(length(int_vars))
        v = check_feasible_flow(all_edges_added, direction, temp_graph, ie, gparams, use_bigm, use_reverse, ndlmo.distmx)
    else
        throw(ArgumentError("Invalid value of use_adaptive $(use_adaptive)"))
    end

    if v == MOI.INFEASIBLE
        #println("Feasibility Tests Failed")
        return Boscia.INFEASIBLE
    end

    #println("Feasibility Tests Passed")
    return Boscia.OPTIMAL
end

function Boscia.is_valid_split(tree::Bonobo.BnBTree, ndlmo::BNDesignLMO, vidx::Int)
    return ndlmo.bounds[vidx, :lessthan] != ndlmo.bounds[vidx, :greaterthan]
end

## Logs
function Boscia.get_BLMO_solve_data(ndlmo::BNDesignLMO)
    return ndlmo.solving_time, 0.0, 0.0
end