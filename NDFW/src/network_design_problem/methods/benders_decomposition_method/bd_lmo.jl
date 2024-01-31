mutable struct BDDesignLMO <: Boscia.BoundedLinearMinimizationOracle
    n::Int64
    graph
    int_vars::Vector{Int}
    bounds::Boscia.IntegerBounds
    gparams::GraphParams
    use_adaptive::Bool
    solving_time::Float64
    distmx::Matrix{Float64}
end

BDDesignLMO(n, graph, int_vars, bounds, gparams, use_adaptive, distmx) = BDDesignLMO(n, graph, int_vars, bounds, gparams, use_adaptive, 0.0, distmx)


function Boscia.compute_extreme_point(
    lmo::BDDesignLMO,
    direction;
    kwargs...
)

    soln = benders_decomposition(lmo, direction)

    return soln
end

function Boscia.build_global_bounds(bdlmo::BDDesignLMO, integer_variables)
    global_bounds = Boscia.IntegerBounds()
    for i in 1:bdlmo.n
        if i in integer_variables
            push!(global_bounds, (i, bdlmo.bounds[i, :lessthan]), :lessthan)
            push!(global_bounds, (i, bdlmo.bounds[i, :greaterthan]), :greaterthan)
        end
    end
    return global_bounds
end

## Read information from problem
function Boscia.get_list_of_variables(bdlmo::BDDesignLMO)
    return bdlmo.n, collect(1:bdlmo.n)
end

# Get list of integer variables, respectively.
function Boscia.get_integer_variables(bdlmo::BDDesignLMO)
    return bdlmo.int_vars
end

function Boscia.get_int_var(bdlmo::BDDesignLMO, cidx)
    return cidx
end

function Boscia.get_lower_bound_list(bdlmo::BDDesignLMO)
    return keys(bdlmo.bounds.lower_bounds)
end

function Boscia.get_upper_bound_list(bdlmo::BDDesignLMO)
    return keys(bdlmo.bounds.upper_bounds)
end

function Boscia.get_bound(bdlmo::BDDesignLMO, c_idx, sense::Symbol)
    @assert sense == :lessthan || sense == :greaterthan
    return bdlmo[c_idx, sense]
end

## Changing the bounds constraints.
function Boscia.set_bound!(bdlmo::BDDesignLMO, c_idx, value, sense::Symbol)
    if sense == :greaterthan
        bdlmo.bounds.lower_bounds[c_idx] = value
    elseif sense == :lessthan
        bdlmo.bounds.upper_bounds[c_idx] = value
    else
        error("Allowed values for sense are :lessthan and :greaterthan.")
    end
end

function Boscia.delete_bounds!(bdlmo::BDDesignLMO, cons_delete)
    # For the cube this shouldn't happen! Otherwise we get unbounded!
    if !isempty(cons_delete)
        error("Trying to delete bounds of the cube!")
    end
end

function Boscia.add_bound_constraint!(bdlmo::BDDesignLMO, key, value, sense::Symbol)
    # Should not be necessary
    return error("Trying to add bound constraints of the cube!")
end

## Checks

function Boscia.is_constraint_on_int_var(bdlmo::BDDesignLMO, c_idx, int_vars)
    return c_idx in int_vars
end

function Boscia.is_bound_in(bdlmo::BDDesignLMO, c_idx, bounds)
    return haskey(bounds, c_idx)
end

function Boscia.is_linear_feasible(bdlmo::BDDesignLMO, v::AbstractVector)
    for i in bdlmo.int_vars
        if !(
            bdlmo.bounds[i, :greaterthan] ≤ v[i] + 1e-6 ||
            !(v[i] - 1e-6 ≤ bdlmo.bounds[i, :lessthan])
        )
            @debug(
                "Vertex entry: $(v[i]) Lower bound: $(bdlmo.bounds[i, :greaterthan]) Upper bound: $(bdlmo.bounds[i, :lessthan]))"
            )
            return false
        end
    end
    return true
end

function Boscia.has_integer_constraint(bdlmo::BDDesignLMO, idx)
    return idx in bdlmo.int_vars
end


###################### Optional
## Safety Functions

function Boscia.build_LMO_correct(bdlmo::BDDesignLMO, node_bounds)
    for key in keys(node_bounds.lower_bounds)
        if !haskey(bdlmo.bounds, (key, :greaterthan)) ||
           bdlmo.bounds[key, :greaterthan] != node_bounds[key, :greaterthan]
            return false
        end
    end
    for key in keys(node_bounds.upper_bounds)
        if !haskey(bdlmo.bounds, (key, :lessthan)) ||
           bdlmo.bounds[key, :lessthan] != node_bounds[key, :lessthan]
            return false
        end
    end
    return true
end

function Boscia.check_feasibility(bdlmo::BDDesignLMO)
    #println("Feasiblity test started")
    for i in bdlmo.int_vars
        if !haskey(bdlmo.bounds, (i, :greaterthan)) || !haskey(bdlmo.bounds, (i, :lessthan))
            return Boscia.UNBOUNDED
        end
        if bdlmo.bounds[i, :greaterthan] > bdlmo.bounds[i, :lessthan]
            return Boscia.INFEASIBLE
        end
    end

    int_vars = bdlmo.int_vars
    gparams = bdlmo.gparams
    ie = length(gparams.init_nodes)
    od_pair_count = gparams.od_pair_count
    total_scenarios = gparams.total_scenarios

    lower_bounds = zeros(length(int_vars))
    upper_bounds = ones(length(int_vars))

    for i in eachindex(int_vars)
        lower_bounds[i] = bdlmo.bounds[int_vars[i], :greaterthan]
        upper_bounds[i] = bdlmo.bounds[int_vars[i], :lessthan]
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

    direction = ones(total_scenarios * (od_pair_count * ie + ie))

    temp_graph = copy(bdlmo.graph)
    if bdlmo.use_adaptive == true
        v = check_feasible_flow(temp_results_y, direction, temp_graph, ie, gparams, bdlmo.distmx)
    elseif bdlmo.use_adaptive == false
        all_edges_added = ones(length(int_vars))
        v = check_feasible_flow(all_edges_added, direction, temp_graph, ie, gparams, bdlmo.distmx)
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

function Boscia.is_valid_split(tree::Bonobo.BnBTree, bdlmo::BDDesignLMO, vidx::Int)
    return bdlmo.bounds[vidx, :lessthan] != bdlmo.bounds[vidx, :greaterthan]
end

## Logs
function Boscia.get_BLMO_solve_data(bdlmo::BDDesignLMO)
    return bdlmo.solving_time, 0.0, 0.0
end

function check_feasible_flow(temp_results_y, direction, graph, ie, gparams, distmx)
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
        
        v = all_or_nothing_stochastic_common(costs, gparams, temp_graph, graph, scen, distmx)
        

        if v == MOI.INFEASIBLE
            return MOI.INFEASIBLE
        end

        orig_v[(first_var+1):(first_var+(od_pair_count*current_edges)+current_edges)] .= v
    end


    final_v = map_to_original(orig_v, temp_graph, gparams, total_scenarios)

    return final_v
end