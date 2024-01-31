using NDFW

ds = "Berlin-Friedrichshain_4"
alg = "CONIC"
frac_rem = 0.01
tn = NDFW.load_ta_network(ds)
total_scenarios = 1

time_limit = 10800
normalize = true
#penalties = [1e1, 1e2, 1e3]
max_fw_iter = 5000
#powers = [1.5]
penalty = 0.0
sres = 1.0
cost_of_expansion = 1000.0
pv = 0
sig_dig = 3
exp_id = 6
rem_cd = true
type_feas_cuts = ["point"]
fraction_removed = frac_rem
seed = Int(round(1000 * (4 + fraction_removed)))
rel_dual_gap = 1e-3
bs_verb = true
fw_verb = true
ls_verb = false
dual_tightening = true
use_bigm = false
use_adaptive = false
use_reverse = false
scale = 0.1
solver = "SCIP"
type_of_nd_constraints = ""

econfig = NDFW.ExperimentConfig(
    ds,
    alg,
    normalize,
    time_limit,
    penalty,
    max_fw_iter,
    pv,
    sig_dig,
    rem_cd,
    type_feas_cuts,
    fraction_removed,
    seed,
    rel_dual_gap,
    bs_verb,
    fw_verb,
    ls_verb,
    dual_tightening,
    total_scenarios,
    use_bigm,
    use_adaptive,
    use_reverse,
    scale,
    type_of_nd_constraints,
    solver
)




soln = NDFW.run_conic_method(cost_of_expansion, sres, econfig)


println("Objective Value: $(soln.result)")
println("Removed Edges Opened: $(soln.new_edges_opened)")
println("Removed edge count: $(soln.new_edge_count)")
println("Solution Vector: $(soln.soln_vector)")

#println("Evaluated solution $(f(soln.soln_vector))")