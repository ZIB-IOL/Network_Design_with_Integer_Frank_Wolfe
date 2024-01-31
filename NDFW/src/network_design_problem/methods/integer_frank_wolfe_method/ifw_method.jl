"""
This function runs the Traffic assignment problem with the Frank Wolfe algorithm 
with and without the network design constraints. When network design constraints 
are present it uses Boscia.jl and when they are not present it uses blended_pairwise_conditional_gradient from the 
FrankWolfe.jl library. Note that the arcs which are to be expanded are determined at random. 
INPUT
ta_data_name : Name of data set from the TransportationNetworks library
enable_network_design : whether to use the network design constraints
base_cost : cost per link of network expansion
res: solution of the tailored frank wolfe algorithm. to be used for scaling. 
OUTPUT:
The above three inputs ans a result which constains the optimal objective function
"""
function network_design(cost_of_expansion, sres, econfig, edges2remove = 0)
    ta_data_name = econfig.ta_data_name
    alg = econfig.alg
    base_cost = cost_of_expansion
    res = sres
    time_limit = econfig.time_limit
    max_fw_iter = econfig.max_fw_iter
    remove_circular_demand = econfig.remove_circular_demand
    type_feas_cuts = econfig.type_feas_cuts
    fraction_removed = econfig.fraction_removed
    bs_verb = econfig.bs_verb
    fw_verb = econfig.fw_verb
    ls_verb = econfig.ls_verb
    dual_tightening = econfig.dual_tightening
    rel_dual_gap = econfig.rel_dual_gap
    total_scenarios = econfig.total_scenarios
    use_bigm = econfig.use_bigm
    use_adaptive = econfig.use_adaptive
    use_reverse = econfig.use_reverse
    scale = econfig.scale
    type_of_nd_constraints = econfig.type_of_nd_constraints
    solver = econfig.solver

    @assert alg == "IFW"

    seed = econfig.seed
    @info "Choosing random seed $(seed)"
    Random.seed!(seed)

    tn, g1, params, gparams = initialize_parameters(ta_data_name, alg, base_cost, remove_circular_demand, total_scenarios, scale)

    #The following lines use the previously defined parameter values to crate the optimization problem. 


    optimizer, con_list = initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios, solver)
    
    cost_of_expansion = []


    # chose 50% of arcs on which allow for network expansion
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
    
    #println("edges to remove: $(new_edge_list)")


    @info "Removed $(new_edge_count) our of $(length(params.edge_list)) edges"
    cost_of_expansion = base_cost * ones(new_edge_count)

    gparams.removed_edges = new_edge_list
    gparams.removed_edges_costs = cost_of_expansion
    

    choose_and_add_network_design_constraints(optimizer, new_edge_list, con_list, g1, params, total_scenarios, type_of_nd_constraints)

    

    lmo = initialize_lmo(alg, optimizer, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)

    if type_of_nd_constraints == "bigM"
        f, grad! = build_f_and_grad(alg, tn, gparams, cost_of_expansion, res, total_scenarios)
    elseif type_of_nd_constraints == "indicator" || type_of_nd_constraints == "both"
        f, grad! = build_f_and_grad_indicator(alg, tn, gparams, cost_of_expansion, res, total_scenarios)
    else
        @error "Invalid type of network design constraints $(type_of_nd_constraints)"
    end

    

    @info "Solving Problem!"

    od_pair_count = gparams.od_pair_count
    no_edges = length(gparams.init_nodes)

    #Choose correct algorithm depending on whether network design constraints are to be used or not. 
    
    x, _, result = Boscia.solve(
        f, 
        grad!, 
        lmo, 
        verbose = bs_verb, 
        time_limit=time_limit, 
        max_fw_iter=max_fw_iter, 
        print_iter = 1,
        fw_verbose = fw_verb,
        dual_tightening=dual_tightening,
        rel_dual_gap=rel_dual_gap,
        #line_search = FrankWolfe.Adaptive(verbose=ls_verb),
        line_search = FrankWolfe.Goldenratio(),
        #line_search=FrankWolfe.AdaptiveZerothOrder(verbose=ls_verb),
        use_postsolve = false
        #strong_convexity = minimum(cost_of_expansion)/2.001
        )

    obj = result[:primal_objective]
    non = result[:number_nodes]
    total_time = result[:total_time_in_sec]
    rel_dual_gap = result[:rel_dual_gap]  
    abs_dual_gap = result[:dual_gap]
    
    total_scenario_vars = od_pair_count * no_edges + no_edges

    x_nd = x[(total_scenario_vars*total_scenarios+1):(total_scenario_vars*total_scenarios+new_edge_count)]
    new_edges_opened = sum(x_nd)

    
    #return [ta_data_name, alg, obj, base_cost, new_edges_opened, new_edge_count, x, new_edge_list, non]
    println("Solution: $(obj)")
    println("new_edge_count: $(new_edge_count)")
    println("new_edges_opened: $(new_edges_opened)")
    println("number of nodes $(non)")
    return Solution(ta_data_name, alg, obj, base_cost, new_edges_opened, 
                new_edge_count, 0, x, new_edge_list, 
                non, total_time, rel_dual_gap, abs_dual_gap)
end


function choose_and_add_network_design_constraints(o, new_edge_list, con_list, g1, params, total_scenarios, type_of_nd_constraints)
    if type_of_nd_constraints == "bigM"
        y = add_network_design_constraints!(
            o,
            new_edge_list,
            con_list,
            g1, 
            params,
            total_scenarios,
            false,
            nothing
        )
    elseif type_of_nd_constraints == "indicator"
        y = add_network_design_constraints_indicator!(
            o,
            new_edge_list,
            con_list,
            g1, 
            params,
            total_scenarios,
            false,
            nothing
        )
    elseif type_of_nd_constraints == "both"
        y = add_network_design_constraints_inverted!(
            o,
            new_edge_list,
            con_list,
            g1, 
            params,
            total_scenarios,
            false, 
            nothing
        )

        add_network_design_constraints_indicator!(
            o,
            new_edge_list,
            con_list,
            g1, 
            params,
            total_scenarios,
            true,
            y
        )
    else
        @error "Invalid type of network design constraints $(type_of_nd_constraints)"
    end
end