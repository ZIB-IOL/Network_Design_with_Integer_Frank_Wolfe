using Test
using Boscia
using NDFW
using LinearAlgebra
using Random
using StatsBase
using JuMP
using FrankWolfe

alg = "IFW"
#ds = "Berlin-Friedrichshain"
#ds = "Berlin-Mitte-Center"
ds = "Braess-Example"
base_cost = 77.2
rem_cd = true
total_scenarios = 1
scale = 0.1
fraction_removed = 0.02
seed = Int(round(1000 * (4 + fraction_removed)))
res = 1.0
penalty = 1000.0
pv = 1.5
use_bigm = false
use_adaptive = false
use_reverse = false
type_feas_cuts = ["point"]
time_limit = 1000
max_fw_iter = 5000
dual_tightening = true
rel_dual_gap_ifw = 1e-4

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

f, grad! = NDFW.build_f_and_grad_indicator(alg, tn, gparams, cost_of_expansion, res, total_scenarios)


type_of_nd_constraints_M = "indicator"

optimizer_M, con_list_M = NDFW.initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios)

NDFW.choose_and_add_network_design_constraints(optimizer_M, new_edge_list, con_list_M, g1, params, total_scenarios, type_of_nd_constraints_M)

MOI.write_to_file(optimizer_M, "test_I.mof.json")

lmo_M = NDFW.initialize_lmo("IFW", optimizer_M, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)


x_ifw_M, _, result_ifw_M = Boscia.solve(
    f,
    grad!,
    lmo_M,
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

obj_ifw_M = result_ifw_M[:primal_objective]

