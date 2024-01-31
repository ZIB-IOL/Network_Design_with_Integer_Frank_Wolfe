using DataFrames
using CSV
using NDFW
using Dates
using Tables
using FrankWolfe
using LinearAlgebra
using Graphs
using Random

ds = "Barcelona"

alg = "SP"
total_scenarios = 1
penalty = 1000.0
pv = 1.5
frac_rem = 0.01
solver = "SCIP"

seed = 4
@info "Choosing random seed $(seed)"
Random.seed!(seed)
base_cost = 0.0
sig_dig = 3
exp_id = 6
rem_cd = true
type_feas_cuts = ["point"]
seed = 1
rel_dual_gap = 1e-4
bs_verb = true
fw_verb = true
ls_verb = false
dual_tightening = true
#total_scenarios = 12
scale = 0.1

tn, g1, params, gparams = NDFW.initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)

println("total_demand: $(sum(gparams.travel_demand[1]))")
#The following lines use the previously defined parameter values to crate the optimization problem. 

optimizer = nothing
lmo = NDFW.initialize_lmo_ta(alg, optimizer, g1.graph, gparams)

f, grad! = NDFW.build_f_and_grad_sp(alg, tn, gparams)
no_edges = length(gparams.init_nodes)
x00 = FrankWolfe.compute_extreme_point(lmo, zeros(no_edges))

storage = zeros(no_edges)
gsp = grad!(storage, x00)

nv = Graphs.nv(g1.graph)
distmx = zeros(nv, nv)

gtn = NDFW.BPR(x00, tn)

@assert norm(gsp - gtn) ≈ 0.0

println("running sp method")
#x01 = FrankWolfe.compute_extreme_point(lmo, gsp)
x01 = NDFW.all_or_nothing_sp(gsp, gparams, g1.graph, 1, gparams.link_dic, distmx)

println("running tn method")
xtn = NDFW.all_or_nothing(gsp, tn, g1.graph, gparams.link_dic)

