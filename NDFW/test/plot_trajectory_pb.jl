using DataFrames
using CSV
using NDFW
using Dates
using Tables
using Random
using StatsBase
using FrankWolfe
using Boscia



ds = "Berlin-Friedrichshain"
#ds = "Anaheim"
#ds = "Berlin-Mitte-Center"
#ds = "Braess-Example"
alg = "BNDLMO-PR"
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
use_reverse = true
alg = "BNDLMO-P"




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


@info "Choosing random seed $(seed)"
Random.seed!(seed)

tn, g1, params, gparams = NDFW.initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)

#The following lines use the previously defined parameter values to crate the optimization problem. 

optimizer, con_list = NDFW.initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios)

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


lmo = NDFW.initialize_lmo(alg, optimizer, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)

f, grad! = NDFW.build_f_and_grad_penalty_method(alg, tn, gparams, cost_of_expansion, res, total_scenarios, use_bigm, penalty, pv)



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
    fw_verbose=fw_verb,
    rel_dual_gap=rel_dual_gap,
    line_search=FrankWolfe.Adaptive(verbose=ls_verb),
    #line_search = FrankWolfe.AdaptiveZerothOrder(verbose=ls_verb),
    dual_tightening=dual_tightening,
    use_postsolve=false
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
    first_var = (scen - 1) * total_scenario_vars
    x_agg = @view(x[(first_var+od_pair_count*no_edges+1):(first_var+od_pair_count*no_edges+no_edges)])

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

bounds_dict = Dict(:list_time => result[:list_time], :list_lb => result[:list_lb], :list_ub => result[:list_ub])

bounds = DataFrame(bounds_dict)

CSV.write("pb_trajectory.csv", Tables.table(bounds), writeheader=false)
