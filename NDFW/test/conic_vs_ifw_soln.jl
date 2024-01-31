using Test
using Boscia
using NDFW
using LinearAlgebra
using Random
using StatsBase
using JuMP
using FrankWolfe

alg = "IFW"
#ds = "Berlin-Friedrichshain_2"
ds = "Berlin-Mitte-Center"
#ds = "Braess-Example"
base_cost = 77.2
rem_cd = true
total_scenarios = 1
scale = 0.1
fraction_removed = 0.02
seed = Int(round(1000 * (4 + fraction_removed)))
res = 386.0
penalty = 1000.0
pv = 1.5
use_bigm = false
use_adaptive = false
use_reverse = false
type_feas_cuts = ["point"]
time_limit = 1000
max_fw_iter = 5000
dual_tightening = true
rel_dual_gap_ifw = 1e-3

@info "Choosing random seed $(seed)"
Random.seed!(seed)

tn, g1, params, gparams = NDFW.initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)

cost_of_expansion = []

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

od_pair_count = gparams.od_pair_count
edge_count = length(gparams.init_nodes)

optimizer, con_list = NDFW.initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios)


o = optimizer

NDFW.add_network_design_constraints!(
    o,
    new_edge_list,
    con_list,
    g1,
    params,
    total_scenarios
)


lmo = NDFW.initialize_lmo("IFW", optimizer, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)

td_2_dest = [maximum([sum(gparams.travel_demand[s][:, z]) for s in 1:total_scenarios]) for z in 1:od_pair_count]
#Total travel demand to each destination

Z_set = collect(1:gparams.od_pair_count)

if gparams.od_pair_count == params.no_nodes
    O_set = collect(1:gparams.od_pair_count)
else
    O_set = collect(gparams.od_pair_count + 1:params.no_nodes)
end


optimizer2 = NDFW.initialize_optimizer()

m, y, x, w = NDFW.build_conic_perspective_model(
    tn,
    cost_of_expansion, 
    params.no_nodes, 
    params.edge_list,
    Z_set,
    O_set, 
    params.node_demands_dict, 
    new_edge_list,
    td_2_dest,
    optimizer2,
    gparams,
    res
)

optimize!(m)

time_conic = solve_time(m)
rel_dual_gap_conic = relative_gap(m)
abs_dual_gap_conic = rel_dual_gap_conic * objective_value(m)

#solution_summary(m)

objval_conic = objective_value(m)
x_val_conic = value.(x)
y_val_conic = value.(y)


f, grad! = NDFW.build_f_and_grad(alg, tn, gparams, cost_of_expansion, res, total_scenarios)

x_ifw, _, result_ifw = Boscia.solve(
    f,
    grad!,
    lmo,
    verbose=true,
    time_limit=time_limit,
    max_fw_iter=max_fw_iter,
    print_iter=1,
    fw_verbose=true,
    dual_tightening=dual_tightening,
    rel_dual_gap=rel_dual_gap_ifw,
    line_search=FrankWolfe.Adaptive(verbose=true),
    #line_search=FrankWolfe.AdaptiveZerothOrder(verbose=ls_verb),
    use_postsolve=false
    #strong_convexity = minimum(cost_of_expansion)/2.001
)

obj_ifw = result_ifw[:primal_objective]
total_time_ifw = result_ifw[:total_time_in_sec]
rel_dual_gap_ifw = result_ifw[:rel_dual_gap]
abs_dual_gap_ifw = result_ifw[:dual_gap]


println("conic objecdive: $(res * objval_conic)")
println("ifw objective: $(res * (obj_ifw - res))")


