"""
Take the current solution. 
Calculate objective 
"""
function optimize_with_penalty(o::Optimizer, direction)
    ie = Graphs.ne(o.graph)
    od_pair_count = o.gparams.od_pair_count
    total_scenarios = o.gparams.total_scenarios

    total_scenario_vars = od_pair_count * ie + ie

    bin_var_costs = direction[(total_scenario_vars*total_scenarios+1):end]
    #res_direction = direction[1:(od_pair_count*ie)+ie]

    lower_bounds, upper_bounds = get_lmo_bounds(o) 
    fixed_vars, relaxed_vars = get_vars_by_state(o)

    # Take the vector of costs and set the positive ones to 0 and negative ones to 1.
    # This is the minimizer. 
    temp_results_y = 1.0 * (min.(0, bin_var_costs) .!= 0)

    int_vars = o.integer_variables


    #We update the computed solution to account for the fixed variables. 
    #Since there are not constraints in the problem all this means is that 
    #we set the fixed variables to their fixed values.
    for i in eachindex(o.gparams.removed_edges)
        if int_vars[i] in fixed_vars
            @assert lower_bounds[i] == upper_bounds[i]

            if upper_bounds[i] == 0.0
                temp_results_y[i] = 0.0
            elseif lower_bounds[i] == 1.0
                temp_results_y[i] = 1.0
            else
                throw(ArgumentError("Incorrect bounds for fixed variable."))
            end
        elseif int_vars[i] in relaxed_vars
            #Do nothing
        else
            throw(ArgumentError("Variable not in fixed or relaxed set."))
        end
    end


    #actually_removed_edges = modify_graph(o, temp_results_y)

    o.curr_binary_soln = temp_results_y
    orig_v = zeros(total_scenario_vars * total_scenarios)

    #new_direction_begin, new_direction_end = obtain_new_directions(o, res_direction, actually_removed_edges, ie)
    for scen = 1:total_scenarios
        first_var = (scen - 1) * total_scenario_vars
        costs = direction[(first_var+1) : (first_var + (od_pair_count * ie)+ie)]


        v, _, _ = all_or_nothing_stochastic_detailed(costs, o.gparams, o.graph, o.true_graph, scen)

        if v == MOI.INFEASIBLE
            o.curr_status = "Infeasible"
            break
        else
            o.curr_status = "Solved"
        end


        
        #orig_v = map_to_original(o, v)
        orig_v[(first_var+1):(first_var+(od_pair_count*ie)+ie)] .= v
    end

    

    #This has to be after mapping v to the original vector
    #restore_graph(o, actually_removed_edges)


    @assert Boscia.is_linear_feasible(o.sopt, [orig_v; temp_results_y])

    # if v == MOI.INFEASIBLE
    #     o.curr_status = "Infeasible"
    # else
    #     o.curr_status = "Solved"
    # end

    term_st = MOI.get(o, MOI.TerminationStatus())

    if term_st ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.SLOW_PROGRESS)
        @error "Unexpected termination: $term_st"
        return o.curr_binary_soln
    end

    o.curr_solution = [orig_v; o.curr_binary_soln]

    @assert Graphs.ne(o.graph) == ie

    return o.curr_solution
end