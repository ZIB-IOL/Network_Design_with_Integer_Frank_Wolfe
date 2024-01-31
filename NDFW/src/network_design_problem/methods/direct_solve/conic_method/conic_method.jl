function run_conic_method(cost_of_expansion, sres, econfig)
    ds = econfig.ta_data_name
    alg = econfig.alg
    base_cost = cost_of_expansion
    res = sres
    time_limit = econfig.time_limit
    max_fw_iter = econfig.max_fw_iter
    rem_cd = econfig.remove_circular_demand
    type_feas_cuts = econfig.type_feas_cuts
    fraction_removed = econfig.fraction_removed
    bs_verb = econfig.bs_verb
    fw_verb = econfig.fw_verb
    ls_verb = econfig.ls_verb
    dual_tightening = econfig.dual_tightening
    rel_dual_gap = econfig.rel_dual_gap
    total_scenarios = econfig.total_scenarios
    use_bigm = econfig.use_bigm
    scale = econfig.scale

    if alg == "CONIC"
        use_pers_form = false
    elseif alg == "CONIC_PERS"
        use_pers_form = true
        alg = "CONIC"
    else
        error("Invalid algorithm")
    end

    tn, g1, params, gparams = initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)

    optimizer = initialize_optimizer()

    seed = econfig.seed
    @info "Choosing random seed $(seed)"
    Random.seed!(seed)

    count_arcs_to_remove = Int(ceil(fraction_removed * length(params.edge_list)))
    new_edge_list = StatsBase.sample(params.edge_list, count_arcs_to_remove, replace=false)
    new_edge_count = length(new_edge_list)
    edge_count = length(params.edge_list)
    


    @info "Removed $(new_edge_count) our of $(length(params.edge_list)) edges"
    cost_of_expansion = base_cost * ones(new_edge_count)

    gparams.removed_edges = new_edge_list
    gparams.removed_edges_costs = cost_of_expansion
    od_pair_count = gparams.od_pair_count
    
    td_2_dest = [maximum([sum(gparams.travel_demand[s][:, z]) for s in 1:total_scenarios]) for z in 1:od_pair_count]
    #Total travel demand to each destination

    Z_set = collect(1:gparams.od_pair_count)
    
    if gparams.od_pair_count == params.no_nodes
        O_set = collect(1:gparams.od_pair_count)
    else
        O_set = collect(gparams.od_pair_count + 1:params.no_nodes)
    end

    if use_pers_form == false
    
        m, y, x, w = build_conic_model(
            tn,
            cost_of_expansion, 
            params.no_nodes, 
            params.edge_list,
            Z_set,
            O_set, 
            params.node_demands_dict, 
            new_edge_list,
            td_2_dest,
            optimizer,
            gparams,
            res
        )

    elseif use_pers_form == true
        m, y, x, w = build_conic_perspective_model(
            tn,
            cost_of_expansion, 
            params.no_nodes, 
            params.edge_list,
            Z_set,
            O_set, 
            params.node_demands_dict, 
            new_edge_list,
            td_2_dest,
            optimizer,
            gparams,
            res
        )
    end


    optimize!(m)

    time_conic = solve_time(m)
    rel_dual_gap = relative_gap(m)
    abs_dual_gap = rel_dual_gap * objective_value(m)

    println(solution_summary(m))
    noc = MOI.get(m, Pajarito.NumberOfCuts())
    println("Number of Cuts: $(MOI.get(m, Pajarito.NumberOfCuts()))")
    println("Time taken: $(time_conic)")
    

    objval = objective_value(m)
    x_val = value.(x)
    y_val = value.(y)

    new_edges_opened_count = Int.(sum(round.(y_val)))

    total_scenario_variables = od_pair_count * edge_count + edge_count
    soln_vector = zeros(total_scenarios * total_scenario_variables + new_edge_count)

    for s in 1:total_scenarios
        first_var = (s - 1) * total_scenario_variables
        vec_x_val = [x_val[e, z, s] for e in params.edge_list for z in Z_set]
        agg_x_val = [sum([x_val[e, z, s] for z in Z_set]) for e in params.edge_list]
        soln_vector[first_var + 1: first_var + od_pair_count * edge_count] = vec_x_val
        soln_vector[first_var+od_pair_count*edge_count+1:first_var+edge_count*(od_pair_count+1)] .= agg_x_val
    end

    vec_y_val = [y_val[e] for e in new_edge_list]

    soln_vector[total_scenario_variables*total_scenarios+1:total_scenario_variables*total_scenarios+new_edge_count] .= vec_y_val

    #selected_edges = [i for i in 1:length(new_edges_opened) if new_edges_opened[i] == 1]
    new_edges_opened = [new_edge_list[i] for i in 1:length(new_edge_list) if vec_y_val[i] == 1]

    non = noc
    
    println(ds)
    println("Algorithm: $(alg)")
    println("Objective value: $(objval)")
    println("Number of nodes: $(non)")
    println("Number of new edges opened: $(new_edges_opened_count)")
    

    return Solution(ds, alg, objval, base_cost, new_edges_opened_count,
        new_edge_count, 0.0, soln_vector, new_edge_list,
        non, time_conic, rel_dual_gap, abs_dual_gap)

    #println("Optimizer")
    #println(optimizer)
    #cuts_identified = backend(m).optimizer.model.optimizer.int_sols_cuts
    #println("Cuts identified: $(cuts_identified)")
    # open("bfmodel.txt", "w") do io
    #     print(io, backend(m).optimizer.model.optimizer.oa_model)
    # end

    return soln
end

function initialize_optimizer()
    optimizer = optimizer_with_attributes(
            Pajarito.Optimizer,
            "oa_solver" => optimizer_with_attributes(
                HiGHS.Optimizer,
                MOI.Silent() => false,
                "mip_feasibility_tolerance" => 1e-8,
                "mip_rel_gap" => 5e-2,
            ),
            "conic_solver" =>
                optimizer_with_attributes(
                    Hypatia.Optimizer, 
                    MOI.Silent() => false,
                    "iter_limit" => 200,
                    "tol_slow" => 1e-5
                    ),
        )
    return optimizer
end