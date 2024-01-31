"""
In this function we use the benders decomposition algorithm along with the shortest path 
problem to solve the network design problem. 
"""
function benders_decomposition(lmo, direction)
    true_graph = copy(lmo.graph)
    gparams = lmo.gparams
    
    ie = Graphs.ne(true_graph)
    nc = Graphs.nv(true_graph)
    od_pair_count = gparams.od_pair_count
    orig_edge_list = [(e.src, e.dst) for e in edges(true_graph)]
    total_scenarios = gparams.total_scenarios
    var_count = total_scenarios * (ie * od_pair_count + ie)
    GRB_ENV = Gurobi.Env()  
    removed_edges = gparams.removed_edges
    int_vars = lmo.int_vars

    distmx = Inf * ones(nc, nc)
    
    # Read current model to identify which variables are free and which are changing.
    lower_bounds = zeros(length(int_vars))
    upper_bounds = ones(length(int_vars))

    for i in eachindex(int_vars)
        lower_bounds[i] = lmo.bounds[int_vars[i], :greaterthan]
        upper_bounds[i] = lmo.bounds[int_vars[i], :lessthan]
    end

    bin_var_costs = direction[int_vars]
    

    res_direction = direction[1:var_count]
    
    # Initialize the master problem
    bopt, yl, eta = initialize_master_model(gparams, bin_var_costs, int_vars, upper_bounds, lower_bounds)

    
    optimality_cuts = []
    feasibility_cuts = []

    max_iter = 1000
    curr_iter = 0.0
    prev_coeff = 0.0
    
    MOI.optimize!(bopt)
    y_vals = MOI.get.(bopt, MOI.VariablePrimal(), yl)
    eta_vals = MOI.get.(bopt, MOI.VariablePrimal(), eta)

    v = 0.0

    #@info "Initial solution: $(y_vals)"
    actually_removed_edges = []
    removed_edge_idxs = []
    orig_v = []

    old_opt_coeff = 0.0

    #Adding back these edges does not change the path
    old_y_vals = copy(y_vals)
    upperbound = Inf
    lowerbound = 0.0
    eps = 1e-4

    normalizing_consts = []
    optimality_cuts = []
    normalizing_const = 0
    original_coeff_list = Dict()

    final_v = zeros(var_count)

    total_scenario_vars = ie * od_pair_count + ie

    status = false
    
    while curr_iter <= max_iter
        #@info "Iteration $(curr_iter) Current Objectiv $(sum(eta_vals))"
        #@info "Removing edges according to solution $(y_vals)."
        #@info "edges in current graph $(Graphs.ne(o.graph))."

        #note that we are not removing edges from the graph. as such the input to the modify graph command 
        #consists of all ones. 
        
        orig_v = zeros(total_scenario_vars * total_scenarios)

        orig_r = Dict()
        orig_t = Dict()
        orig_s = Dict()


        #@info "edges in current graph $(Graphs.ne(o.graph))."

        new_direction = res_direction

        all_results = []
        for scen in 1:total_scenarios
            first_var = (scen - 1) * (ie * od_pair_count + ie)

            result_vals, new_optimality_cuts, new_feasibility_cuts, vr, vt, vs = single_scenario_cuts(scen, true_graph, orig_edge_list, gparams, new_direction, yl, eta[scen], y_vals, GRB_ENV)
            optimality_cuts = [optimality_cuts; new_optimality_cuts]
            feasibility_cuts = [feasibility_cuts; new_feasibility_cuts]
            push!(all_results, result_vals)

            if result_vals != MOI.INFEASIBLE
                orig_r[scen] = vr
                orig_t[scen] = vt
                orig_s[scen] = vs
            end
        end

        # TODO
        # check for optimality
        if any(all_results .== MOI.INFEASIBLE) == false
            status = check_for_optimality(lowerbound, gparams, orig_r, orig_t, orig_s, new_direction, bin_var_costs, y_vals, eps)
        end


        if curr_iter >= max_iter || status == true
            #@info result_vals
            #@info "Reached max iterations. Returning current solution."
            #@info "Current solution: $(y_vals)"
            #@info "Current objective: $(dot(y_vals, bin_var_costs) + sum(eta_vals))"
            #@info MOI.get.(bopt, MOI.VariablePrimal(), yl)
            #@info MOI.get(bopt, MOI.TerminationStatus())

            temp_graph = copy(true_graph)
            modified_graph, curr_removed_edges = modify_graph(temp_graph, y_vals, gparams.removed_edges)
            modified_costs = modify_cost_vector(new_direction, gparams, curr_removed_edges, ie)
            new_edge_list = [(edge.src, edge.dst) for edge in edges(modified_graph)]
            mie = Graphs.ne(modified_graph)
            
            for scen in 1:total_scenarios
                curr_first_var = (scen - 1) * (mie * od_pair_count + mie)
                orig_first_var = (scen - 1) * (ie * od_pair_count + ie)
                costs = modified_costs[curr_first_var+1:(curr_first_var+mie*od_pair_count+mie)]
                
                v, _, _, _, _ = all_or_nothing_stochastic_detailed(costs, gparams, modified_graph, true_graph, scen, distmx)

                fv = map_to_original(v, temp_graph, gparams, 1)

                orig_v[(orig_first_var+1):(orig_first_var+(od_pair_count*ie)+ie)] .= fv

            end

            #final_v = map_to_original(orig_v, temp_graph, gparams, total_scenarios)

            final_v = orig_v

            obj_val = 0.0
            for scen in 1:total_scenarios
                first_var = (scen - 1) * total_scenario_vars
                odpie = od_pair_count * ie
                obj_val += dot(new_direction[first_var + odpie + 1:first_var + odpie + ie], final_v[first_var + odpie + 1:first_var + odpie + ie])
            end

            break

        else
            # add cuts to problem

            for c in optimality_cuts
                MOI.add_constraint(bopt, c[1][1], c[1][2])
            end

            for c in feasibility_cuts
                MOI.add_constraint(bopt, c[1], c[2])
            end

            #MOI.write_to_file(bopt, "master2.lp")
            MOI.optimize!(bopt)
            @info MOI.get(bopt, MOI.TerminationStatus())
            
            eta_vals = MOI.get.(bopt, MOI.VariablePrimal(), eta)

            y_vals = MOI.get.(bopt, MOI.VariablePrimal(), yl)
            
            #curr_normalizing_const = length(normalizing_consts) == 0 ? normalizing_const : minimum(normalizing_consts)
            #curr_normalizing_const = 0.0
            lowerbound = sum([bin_var_costs[i] * y_vals[i] for i in eachindex(y_vals)]) + sum(eta_vals)
            #@assert norm(y_vals - round.(y_vals)) < 1e-10  
            y_vals = round.(y_vals)

            if any(all_results .== MOI.INFEASIBLE)
                edges_added = findall(x->x==1.0, y_vals)
                #@info "edges added: $(gparams.removed_edges[edges_added])"
            end

            curr_iter += 1
            #@info "edges in new solution $(sum(y_vals))"
        end

        #restore_graph(o, actually_removed_edges)
    end
    
    #@info "Used $(curr_iter) iterations."

    #@info "length of final solution is $(length(v))"
    #@info "Final solution: $(MOI.get.(bopt, MOI.VariablePrimal(), yl))" 
    #@info "Final objective: $(eta_val)"

    @assert all([final_v[edge_count*gparams.od_pair_count+re] == 0.0 for re in removed_edge_idxs])
    @assert length(final_v) == total_scenarios * (Graphs.ne(true_graph) * gparams.od_pair_count + Graphs.ne(true_graph))

    zero_edges_zero = []
    for re in actually_removed_edges
        edge_idx = gparams.link_dic[re[1], re[2]]
        edge_val = final_v[ie*gparams.od_pair_count+edge_idx]
        push!(zero_edges_zero, edge_val == 0.0)
    end
    @assert all(zero_edges_zero)

    
    curr_binary_soln = y_vals

    removed_edge_idxs = [gparams.link_dic[eg[1],eg[2]] for eg in actually_removed_edges]

    # if v == MOI.INFEASIBLE
    #     o.curr_status = "Infeasible"
    # else
    #     o.curr_status = "Solved"
    # end

    # term_st = MOI.get(o, MOI.TerminationStatus())

    # if term_st ∉ (MOI.OPTIMAL, MOI.ALMOST_OPTIMAL, MOI.SLOW_PROGRESS)
    #     @error "Unexpected termination: $term_st"
    #     return o.curr_binary_soln
    # end

    curr_solution = [final_v; curr_binary_soln]

    obj_val = 0.0
    for scen in 1:total_scenarios
        first_var = (scen - 1) * total_scenario_vars
        odpie = od_pair_count * ie
        obj_val += dot(direction[first_var + odpie + 1:first_var + odpie + ie], final_v[first_var + odpie + 1:first_var + odpie + ie])
    end

    @assert Graphs.ne(true_graph) == ie

    return curr_solution
end

function single_scenario_cuts(scen, modified_graph, edge_list, gparams, new_direction, yl, eta, y_vals, GRB_ENV)
    od_pair_count = gparams.od_pair_count
    edge_count = Graphs.ne(modified_graph)
    node_count = Graphs.nv(modified_graph)

    first_var = (scen - 1) * (edge_count * od_pair_count + edge_count)

    res_direction = new_direction[first_var + 1:(first_var + edge_count * od_pair_count + edge_count)]

    result_vals, vr, vt, vs = all_or_nothing_LP(scen, res_direction, gparams, y_vals, node_count, edge_list, GRB_ENV)

    feasibility_cuts = []
    optimality_cuts = []
            
    if result_vals == MOI.INFEASIBLE 
        #result_vals2 = all_or_nothing_detailed(costs, o.gparams, o.graph, o.true_graph)

        #bopt = add_feasibility_cuts(o, bopt, yl, state, source, costs, orig_edge_list)
        if all(vs .== 0)
            @info "Infeasible solution. Adding bridge feasibility cuts."
            status, state, source = all_or_nothing_detailed(costs, gparams, graph, true_graph)
            @info "network sol status $(status)"
            new_feasibility_cuts = add_feasibility_cuts(gparams, yl, y_vals, state, source, costs, orig_edge_list, type_feas_cuts)
        else
            @info "Infeasible solution. Adding subproblem feasibility cuts."
            new_feasibility_cuts = add_feasibility_cuts(gparams, yl, vr, vt, vs, scen)
        end

        feasibility_cuts = [feasibility_cuts; new_feasibility_cuts]
        #bopt = add_feasibility_cuts(o, bopt, yl, vr, vt, vs)
    else        
        #r, t, s = result_vals   
        r = vr
        t = vt
        s = vs

        # check for optimality
        # check_for_optimality()

        c = add_optimality_cuts(scen, gparams, yl, eta, r, s, t)

        # Normalize the cut by reducing the RHS
        # normalize_cut(c, normalizing_consts)

        push!(optimality_cuts, c)

    end

    return result_vals, optimality_cuts, feasibility_cuts, vr, vt, vs
end

"""
This function initializes the model of the master problem. 
INPUTS:
    o: Optimizer object
    bin_var_costs: cost of each binary variable
    fixed_vars: integer variables that are fixed to 0 or 1
    relaxed_vars: integer variables that are relaxed
    upper_bounds: upper bounds of the integer variables
    lower_bounds: lower bounds of the integer variables
"""
function initialize_master_model(gparams, bin_var_costs, int_vars, upper_bounds, lower_bounds)
    bopt = SCIP.Optimizer()
    MOI.set(bopt, MOI.Silent(), true)
    var_count = length(gparams.removed_edges)
    total_scenarios = gparams.total_scenarios

    yl = MOI.add_variables(bopt, var_count)
    eta = MOI.add_variables(bopt, total_scenarios)

    int_vars = int_vars

    for i in eachindex(int_vars)
        MOI.add_constraint(bopt, yl[i], MOI.ZeroOne())

        if lower_bounds[i] == upper_bounds[i]

            if upper_bounds[i] == 0.0
                MOI.add_constraint(bopt, yl[i], MOI.LessThan(0.0))
                MOI.add_constraint(bopt, yl[i], MOI.GreaterThan(0.0))
            elseif lower_bounds[i] == 1.0
                MOI.add_constraint(bopt, yl[i], MOI.LessThan(1.0))
                MOI.add_constraint(bopt, yl[i], MOI.GreaterThan(1.0))
            else
                throw(ArgumentError("Incorrect bounds for fixed variable."))
            end

        elseif lower_bounds[i] < upper_bounds[i]
            MOI.add_constraint(bopt, yl[i], MOI.LessThan(1.0))
            MOI.add_constraint(bopt, yl[i], MOI.GreaterThan(0.0))

        else
            throw(ArgumentError("Variable not in fixed or relaxed set."))
        end

        MOI.set(bopt, MOI.VariableName(), yl[i], "y_$(i)")
    end

    for i in eachindex(eta)
        MOI.set(bopt, MOI.VariableName(), eta[i], "eta_$(i)")
        MOI.add_constraint(bopt, eta[i], MOI.GreaterThan(-sum(bin_var_costs)))
    end


    MOI.set(bopt, MOI.ObjectiveSense(), MOI.MIN_SENSE)

    design_term = MOI.ScalarAffineTerm.(bin_var_costs, yl)
    stoch_term = MOI.ScalarAffineTerm.(ones(total_scenarios), eta)
    objective_terms = [design_term; stoch_term]

    MOI.set(
        bopt,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(objective_terms, 0.0)
    )

    return bopt, yl, eta
end

"""
we want to check if we are running and optimization problem or a feasibility test
The derivatives of the initial variables must be 0 for the optimality problem
"""
function check_if_feasibility_problem(o, direction)
    ie = Graphs.ne(o.graph)
    od_pair_count = o.gparams.od_pair_count

    if all(direction[1:(od_pair_count*ie)] .> 0.0)
        println("Running feasibility test")
        feasibility_test_flag = true
    else
        println("Solving optimization problem")
        feasibility_test_flag = false
    end
    return feasibility_test_flag
end

