using Test
using Boscia
using NDFW
using ForwardDiff
using LinearAlgebra
using Random
using StatsBase
using DelimitedFiles
import MathOptInterface as MOI

alg = "BNDLMO-P"
ds = "Anaheim"
base_cost = 1407.0544087491
rem_cd = true
total_scenarios = 5
scale = 0.1
fraction_removed = 0.02
seed = Int(round(1000 * (4 + fraction_removed)))
res = 1286047.73/100
penalty = 1000.0
pv = 1.5
type_feas_cuts = ["point"]
use_bigm = false

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

optimizer, con_list = NDFW.initialize_model("IFW", tn, g1, params, gparams, type_feas_cuts, total_scenarios)


o = optimizer

NDFW.add_network_design_constraints!(
    o,
    new_edge_list,
    con_list,
    g1,
    params,
    total_scenarios
)


lmo = NDFW.initialize_lmo("IFW", optimizer, g1.graph, gparams, use_bigm)

blmo = NDFW.initialize_lmo("BNDLMO-P", optimizer, g1.graph, gparams, use_bigm)

MOI.write_to_file(lmo.o, "anaheim_model2.lp")


var_count = total_scenarios * (od_pair_count * edge_count + edge_count) + new_edge_count
total_scenario_variables = od_pair_count * edge_count + edge_count

# x = rand(var_count)
# for scen in 1:total_scenarios
#     first_var = (scen - 1) * total_scenario_variables

#     for i in 1:edge_count
#         curr_x_flows = [x[first_var + (z-1) * edge_count + i] for z in 1:od_pair_count]

#         x[first_var + od_pair_count * edge_count + i] = sum(curr_x_flows)
#     end
# end

#x[var_count-new_edge_count+1:var_count] = rand([0, 1], new_edge_count)
#x[var_count-new_edge_count+1:var_count] = ones(new_edge_count)
xi = readdlm("./examples/blmo_test_v4/solution_Anaheim_missing_missing_IFW_0.02_5.csv", ',', Float64)
xb = readdlm("/Users/shakart/Downloads/solution_Anaheim_1000.0_1.5_BNDLMO-P_0.02_5.csv", ',', Float64)

vv = readdlm("./examples/vector_v.csv", ',', Float64)
vx = readdlm("./examples/vector_x.csv", ',', Float64)
vg = readdlm("./examples/vector_gradient.csv", ',', Float64)


alg = "IFW"
fifw, gifw =  NDFW.build_f_and_grad(alg, tn, gparams, cost_of_expansion, res, total_scenarios)

alg = "BNDLMO-P"
fblmo, gblmo = NDFW.build_f_and_grad_penalty_method(alg, tn, gparams, cost_of_expansion, res, total_scenarios, use_bigm, penalty, pv)

#auto_grad = ForwardDiff.gradient(f, x)


#grad = grad!(storage, x)
fi = fifw(xi)
fb = fblmo(xb)

storage = zeros(var_count)
gifw(storage, xi)
gi = copy(storage)

storage = zeros(var_count)
gblmo(storage, xb)
gb = copy(storage)

# vb = Boscia.compute_extreme_point(blmo, gb)
vi = Boscia.compute_extreme_point(lmo, gi)


# #@assert fi ≈ fb
# ubl = [MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex,MathOptInterface.LessThan{Float64}}(i) for i in 178231:178249]
# lbl = [MathOptInterface.ConstraintIndex{MathOptInterface.VariableIndex,MathOptInterface.GreaterThan{Float64}}(i) for i in 178231:178249]

# ubnds = [Boscia.get_bound(lmo, u, :lessthan).upper for u in ubl]
# lbnds = [Boscia.get_bound(lmo, l, :greaterthan).lower for l in lbl]