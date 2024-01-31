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
function traffic_assignment_lp(econfig)
    ta_data_name = econfig.ta_data_name
    alg = econfig.alg
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

    @assert alg == "LP"
    @assert total_scenarios == 1

    seed = econfig.seed
    @info "Choosing random seed $(seed)"
    Random.seed!(seed)
    base_cost = 0.0

    tn, g1, params, gparams = initialize_parameters(ta_data_name, alg, base_cost, remove_circular_demand, total_scenarios, scale)

    #The following lines use the previously defined parameter values to crate the optimization problem. 


    optimizer, con_list = initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios, solver)
    

    lmo = initialize_lmo_ta(alg, optimizer, g1.graph, gparams)

    f, grad! = build_f_and_grad_TA(alg, tn, gparams, total_scenarios)
    
    

    @info "Solving Problem!"

    od_pair_count = gparams.od_pair_count
    no_edges = length(gparams.init_nodes)

    total_vars = od_pair_count * no_edges + no_edges


    x00 = FrankWolfe.compute_extreme_point(lmo, zeros(total_vars))

    #Choose correct algorithm depending on whether network design constraints are to be used or not. 
    
    total_time = @elapsed result = FrankWolfe.blended_pairwise_conditional_gradient(
        f, 
        grad!, 
        lmo, 
        x00,
        verbose = bs_verb, 
        max_iteration=max_fw_iter, 
        print_iter = 1,
        line_search = FrankWolfe.Adaptive(verbose=ls_verb),
        epsilon=1e-3 * sum(tn.travel_demand),
        )
    
    obj = result.primal
    abs_dual_gap = result.dual_gap
    rel_dual_gap = abs_dual_gap/obj
    x = result.x

    #return [ta_data_name, alg, obj, base_cost, new_edges_opened, new_edge_count, x, new_edge_list, non]
    println("Solution: $(obj)")
    return Solution(ta_data_name, alg, obj, 0.0, 0.0, 
                0.0, 0, x, [], 
                0, total_time, rel_dual_gap, abs_dual_gap)
end


