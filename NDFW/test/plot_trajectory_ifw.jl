using Test
using Boscia
using NDFW
using LinearAlgebra
using Random
using StatsBase
using JuMP
using FrankWolfe
using CSV
using Tables
using DataFrames

ds = "Berlin-Friedrichshain"
#ds = "Anaheim"
#ds = "Berlin-Mitte-Center"
#ds = "Braess-Example"
#alg = "BNDLMO-PR"
alg = "IFW"
total_scenarios = 1
penalty = 1000.0
pv = 1.5
frac_rem = 0.01
solver="SCIP"
base_cost = 1181.72069
res = 1.0

use_adaptive = true
type_of_nd_constraints = ""
filealg = alg


use_bigm = false
use_reverse = false
type_of_nd_constraints = "bigM"



println("Running dataset $(ds)")

time_limit = 10800
normalize = true
#penalties = [1e1, 1e2, 1e3]
max_fw_iter = 5000
#powers = [1.5]

sig_dig = 3
exp_id = 6
rem_cd = true
type_feas_cuts = ["point"]
fraction_removed = frac_rem
seed = Int(round(1000 * (4 + fraction_removed)))
rel_dual_gap = 5e-2
bs_verb = true
fw_verb = true
ls_verb = false
dual_tightening = true
#total_scenarios = 12
scale = 0.1
ta_data_name = ds
remove_circular_demand = rem_cd
@info "Choosing random seed $(seed)"
Random.seed!(seed)

tn, g1, params, gparams = NDFW.initialize_parameters(ta_data_name, alg, base_cost, remove_circular_demand, total_scenarios, scale)

#The following lines use the previously defined parameter values to crate the optimization problem. 


optimizer, con_list = NDFW.initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios, solver)

cost_of_expansion = []


# chose 50% of arcs on which allow for network expansion
@info "Choosing random seed $(seed)"
Random.seed!(seed)

#arcs_to_be_expanded = [rand() > 0.5 for i in 1:length(edge_list)]
count_arcs_to_remove = Int(ceil(fraction_removed * length(params.edge_list)))
new_edge_list = StatsBase.sample(params.edge_list, count_arcs_to_remove, replace=false)
new_edge_count = length(new_edge_list)


#println("edges to remove: $(new_edge_list)")


@info "Removed $(new_edge_count) our of $(length(params.edge_list)) edges"
cost_of_expansion = base_cost * ones(new_edge_count)

gparams.removed_edges = new_edge_list
gparams.removed_edges_costs = cost_of_expansion


NDFW.choose_and_add_network_design_constraints(optimizer, new_edge_list, con_list, g1, params, total_scenarios, type_of_nd_constraints)



lmo = NDFW.initialize_lmo(alg, optimizer, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)

if type_of_nd_constraints == "bigM"
    f, grad! = NDFW.build_f_and_grad(alg, tn, gparams, cost_of_expansion, res, total_scenarios)
elseif type_of_nd_constraints == "indicator" || type_of_nd_constraints == "both"
    f, grad! = NDFW.build_f_and_grad_indicator(alg, tn, gparams, cost_of_expansion, res, total_scenarios)
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
    verbose=bs_verb,
    time_limit=time_limit,
    max_fw_iter=max_fw_iter,
    print_iter=1,
    fw_verbose=fw_verb,
    dual_tightening=dual_tightening,
    rel_dual_gap=rel_dual_gap,
    line_search=FrankWolfe.Adaptive(verbose=ls_verb),
    #line_search=FrankWolfe.AdaptiveZerothOrder(verbose=ls_verb),
    use_postsolve=false
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

bounds_dict = Dict(:list_time => result[:list_time], :list_lb => result[:list_lb], :list_ub => result[:list_ub])

bounds = DataFrame(bounds_dict)

CSV.write("ifw_trajectory.csv", bounds)
