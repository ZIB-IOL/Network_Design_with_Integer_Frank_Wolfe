function run_bounded_lmo_penalty_approach(cost_of_expansion, sres, econfig)

    alg = econfig.alg
    ds = econfig.ta_data_name
    base_cost = cost_of_expansion
    res = sres
    time_limit = econfig.time_limit
    penalty = econfig.time_limit
    max_fw_iter = econfig.max_fw_iter
    pv = econfig.pv
    rem_cd = econfig.remove_circular_demand
    type_feas_cuts = econfig.type_feas_cuts
    fraction_removed = econfig.fraction_removed
    rel_dual_gap = econfig.rel_dual_gap
    bs_verb = econfig.bs_verb
    fw_verb = econfig.fw_verb
    ls_verb = econfig.ls_verb
    dual_tightening = econfig.dual_tightening
    total_scenarios = econfig.total_scenarios
    use_bigm = econfig.use_bigm
    use_adaptive = econfig.use_adaptive
    use_reverse = econfig.use_reverse
    scale = econfig.scale

    seed = econfig.seed
    @info "Choosing random seed $(seed)"
    Random.seed!(seed)

    tn, g1, params, gparams = initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)

    #The following lines use the previously defined parameter values to crate the optimization problem. 

    optimizer, con_list = initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios)

    cost_of_expansion = []

    seed = econfig.seed
    @info "Choosing random seed $(seed)"
    Random.seed!(seed)
    
    #arcs_to_be_expanded = [rand() > 0.5 for i in 1:length(edge_list)]
    count_arcs_to_remove = Int(ceil(fraction_removed * length(params.edge_list)))
    new_edge_list = StatsBase.sample(params.edge_list, count_arcs_to_remove, replace=false)
    new_edge_count = length(new_edge_list)
        

    @info "Removed $(new_edge_count) our of $(length(params.edge_list)) edges"
    cost_of_expansion = base_cost * ones(new_edge_count)

    gparams.removed_edges = new_edge_list
    gparams.removed_edges_costs = cost_of_expansion


    lmo = initialize_lmo(alg, optimizer, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)

    f, grad! = build_f_and_grad_penalty_method(alg, tn, gparams, cost_of_expansion, res, total_scenarios, use_bigm, penalty, pv)

        

    @info "Solving Problem!"

    od_pair_count = gparams.od_pair_count
    no_edges = length(gparams.init_nodes)

    #Choose correct algorithm depending on whether network design constraints are to be used or not. 
    x, _, result = Boscia.solve(
        f, 
        grad!, 
        lmo, 
        verbose=bs_verb,
        time_limit=time_limit,
        max_fw_iter=max_fw_iter,
        fw_verbose = fw_verb,
        rel_dual_gap=rel_dual_gap,
        #line_search = FrankWolfe.Adaptive(verbose=ls_verb),
        line_search = FrankWolfe.Goldenratio(),
	    #line_search = FrankWolfe.AdaptiveZerothOrder(verbose=ls_verb),
        dual_tightening = dual_tightening,
        use_postsolve = false
        #strong_convexity = minimum(cost_of_expansion)/2.001
        )

    @info "Solved Problem!"
    
    # for k in keys(result)
    #     @info "$k: $(result[k])"
    # end

    # for v in lmo.int_vars
    #    println("$(lmo.bounds[v,:greaterthan]) <= x[$(v)] <= $(lmo.bounds[v,:lessthan])")
    # end
        
    obj = result[:primal_objective]
    non = result[:number_nodes]
    total_time = result[:total_time_in_sec]
    abs_dual_gap = result[:dual_gap]
    rel_dual_gap = result[:rel_dual_gap]

    total_scenario_vars = od_pair_count * no_edges + no_edges

    x_nd = @view(x[(total_scenario_vars*total_scenarios+1):(total_scenario_vars*total_scenarios+new_edge_count)])

    bigM = tn.total_od_flow
    removed_edge_idx = [gparams.link_dic[removed_edge[1], removed_edge[2]] for removed_edge in gparams.removed_edges]

    violation = []
    for scen = 1:total_scenarios
        first_var = (scen-1)*total_scenario_vars
        x_agg = @view(x[(first_var + od_pair_count*no_edges+1):(first_var + od_pair_count*no_edges+no_edges)])
        
        # for j in eachindex(x_nd)
        #     cons_violation = max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0)
        # end
        
        cons_violation = [max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0) for j in eachindex(x_nd)]
        #println("Violation $(cons_violation)")
        append!(violation, maximum(cons_violation))
    end

    max_violation = maximum(violation)


    #@info "Solved Problem! Objective value: $(result)"
    @info "Max Constraint Violation is $(max_violation)"
    #@info "Constraint Wise Violation: $(cons_violation)"
    @info "Optimal Solution is $(x_nd)"

    return Solution(ds, alg, obj, base_cost, sum(x_nd), 
                    new_edge_count, max_violation, x, new_edge_list,
                    non, total_time, rel_dual_gap, abs_dual_gap)
end
