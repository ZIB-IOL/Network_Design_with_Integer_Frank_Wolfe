struct ExprEvaluator <: MOI.AbstractNLPEvaluator
    constraints::Vector{Expr}
end

function MOI.initialize(d::ExprEvaluator, requested_features::Vector{Symbol})
    if requested_features != [:ExprGraph]
        error("Only supports expression graph!")
    end
end

MOI.constraint_expr(evaluator::ExprEvaluator, i) = evaluator.constraints[i]

"""
This function solves the entire minimizaiton problem using SCIP
"""
function direct_scip_solve_stochastic(econfig, cost_of_expansion, sres, edges2remove=0)
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
    type_of_nd_constraints = econfig.type_of_nd_constraints
    solver = econfig.solver
    
    sf = 1 # scaling factor
    
    base_cost = base_cost / sf

    lower_bound = -Inf
    upper_bound = 0.0

    type_feas_cuts = ["bridge"]


    seed = econfig.seed
    @info "Choosing random seed $(seed)"
    Random.seed!(seed)


    tn, g1, params, gparams = initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)

    #The following lines use the previously defined parameter values to crate the optimization problem. 

    optimizer, con_list = initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios, solver)

    cost_of_expansion = []

    seed = econfig.seed
    @info "Choosing random seed $(seed)"
    Random.seed!(seed)

    if edges2remove == 0
        #arcs_to_be_expanded = [rand() > 0.5 for i in 1:length(edge_list)]
        count_arcs_to_remove = Int(ceil(fraction_removed * length(params.edge_list)))
        new_edge_list = StatsBase.sample(params.edge_list, count_arcs_to_remove, replace=false)
        new_edge_count = length(new_edge_list)
    else
        new_edge_list = edges2remove
        new_edge_count = length(new_edge_list)
    end


    @info "Removed $(new_edge_count) our of $(length(params.edge_list)) edges"
    cost_of_expansion = base_cost * ones(new_edge_count)

    gparams.removed_edges = new_edge_list
    gparams.removed_edges_costs = cost_of_expansion


    choose_and_add_network_design_constraints(optimizer, new_edge_list, con_list, g1, params, total_scenarios, type_of_nd_constraints)


    od_pair_count = gparams.od_pair_count
    edge_count = length(gparams.init_nodes)
    new_edge_count = length(gparams.removed_edges)
    var_cnt = edge_count * od_pair_count
    total_scenario_vars = od_pair_count*edge_count + edge_count

    x = Dict()
    x_agg = Dict()

    for scen in 1:total_scenarios
        first_var = (scen - 1) * total_scenario_vars
        x[scen] = [MOI.VariableIndex(i) for i = (first_var + 1):(first_var + var_cnt)]
        x_agg[scen] = [MOI.VariableIndex(i) for i = (first_var+var_cnt+1):(first_var+var_cnt+edge_count)]
    end

    x_nd = [MOI.VariableIndex(i) for i = (total_scenario_vars*total_scenarios+1):(total_scenario_vars*total_scenarios+new_edge_count)]
    
    
    
    tb = MOI.add_variables(optimizer, edge_count * total_scenarios)
    obj_var = MOI.add_variable(optimizer)

    println("Selected variables")

    #x_agg = @view(x[(od_pair_count*edge_count+1):(od_pair_count*edge_count+edge_count)])

    terms = Dict()

    for scen in 1:total_scenarios
        terms[scen] = Expr[]
        for i in eachindex(x_agg[scen])
            const_term = tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
            om = tn.free_flow_time[i]
            ia = tn.b[i]
            powt = tn.power[i] + 1
            denom = (tn.capacity[i]^tn.power[i]) * (tn.power[i] + 1)

            term = :($(const_term) 
                    + $(om) * (x[$(x_agg[scen][i])] + $(ia) * (x[$(x_agg[scen][i])]^($powt)) / $(denom))
                        - x[$(tb[i])] <= 0.0
                    )


            append!(terms[scen], [term])
        end
    end


    nlpb = MOI.NLPBlockData(
        [MOI.NLPBoundsPair(lower_bound, upper_bound) for i in 1:params.no_edges],
        ExprEvaluator(terms),
        false
    )

    MOI.set(optimizer, MOI.NLPBlock(), nlpb)

    println("Initialized non linear block")

    MOI.add_constraint(optimizer,
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([1.0 * ones(params.no_edges); -1.0 * sf], [tb; obj_var]),
            0.0,
        ),
        MOI.LessThan(0.0),
    )


    # if alg == "IFW" || alg == "NLMO-IFW"
    #     x_nd = @view(x[((od_pair_count+1)*edge_count+1):((od_pair_count+1)*edge_count+new_edge_count)])
    #     for i in eachindex(x_nd)
    #         sum += cost_of_expansion[i] * x_nd[i]
    #     end
    # end

    MOI.set(
        optimizer,
        MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
        MOI.ScalarAffineFunction(
            MOI.ScalarAffineTerm.([1.0; cost_of_expansion], [obj_var; x_nd]),
            0.0,
        ),
    )

    println("Initialized objective function")

    MOI.set(optimizer, MOI.ObjectiveSense(), MOI.MIN_SENSE)
    MOI.set(optimizer, MOI.Silent(), false)

    #MOI.write_to_file(optimizer, "$(ds)_model.mps")
    #SCIP.LibSCIP.SCIPwriteLP(optimizer, "model.cip")

    println("starting optimization")

    total_time = @elapsed MOI.optimize!(optimizer)

    objective_value = MOI.get(optimizer, MOI.ObjectiveValue())
    new_arcs = MOI.get(optimizer, MOI.VariablePrimal(), x_nd)
    arc_flows = MOI.get(optimizer, MOI.VariablePrimal(), x_agg)
    new_edges_opened = sum(new_arcs)

    println(arc_flows)
    println(new_arcs)
    println(new_edge_list)

    return Solution(ds, alg, objective_value, base_cost, new_edges_opened,
                    new_edge_count, 0, [], new_edge_list,
                    0, total_time, 0, 0)

    return [ds, alg, objective_value, base_cost, new_edges_opened, new_edge_count]

end
