using DataFrames
using CSV
using NDFW
using Dates
using Tables



#ds = "Berlin-Friedrichshain"
ds = "Anaheim"
#ds = "Berlin-Mitte-Center"
#ds = "Berlin-Tiergarten"
#ds = "Berlin-Prenzlauerberg-Center"
#ds = "Braess-Example"
#alg = "BNDLMO-PR"
alg = "TA"
total_scenarios = 1
penalty = 1000.0
pv = 1.5
frac_rem = 0.01
solver = "SCIP"


use_adaptive = true
type_of_nd_constraints = ""
filealg = alg

use_bigm = false
use_reverse = false




println("Running dataset $(ds)")

time_limit = 10800
normalize = true
#penalties = [1e1, 1e2, 1e3]
max_fw_iter = 250
#powers = [1.5]

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
#total_scenarios = 12
scale = 0.1



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


results = DataFrame(Name=String[],
    Nodes=Int64[],
    Edges=Int64[],
    Algorithm=String[],
    Result=Float64[],
    BaseCost=Any[],
    Penalty=Float64[],
    PowerCoeff=Float64[],
    FractionRemoved=Float64[],
    RemovedEdgesReopened=Any[],
    EdgesRemoved=Any[],
    TotalConstraintViolation=Float64[],
    Time=Float64[],
    NumberOfNodes=Int64[],
    TimePerNode=Float64[],
    RelDualGap=Float64[],
    AbsDualGap=Float64[],
    TotalScenarios=Int64[],
    Solver=String[]
)

allowmissing!(results)

date_today = Date(string(today()), dateformat"y-m-d")

folder_name = "ta_test"
#folder_name = "experiment_$(Dates.format(date_today, "ddmmyy"))_$(exp_id)" 

isdir(folder_name) || mkdir(folder_name)

results_file = "$(folder_name)/results_$(ds)_$(penalty)_$(pv)_$(filealg)_$(fraction_removed)_$(total_scenarios)_$(solver).csv"

@info "Evaluating data set $(ds)"
tn = NDFW.load_ta_network(ds)

@info "Running tailored FW algorithm"
t = @elapsed _, _, res = NDFW.ta_frank_wolfe(tn)
#alg = "TN"
@info "Obtained solution:" round(res, digits=sig_dig)
output_row = [ds,
    tn.number_of_nodes,
    tn.number_of_links,
    "TN",
    round(res, digits=sig_dig),
    missing,
    missing,
    missing,
    missing,
    missing,
    missing,
    missing,
    round(t, digits=sig_dig),
    missing,
    missing,
    missing,
    missing,
    missing,
    missing
]

push!(results, output_row)
CSV.write(results_file, results)

#for alg in ["FW", "NLMO-FW", "IFW", "NLMO-IFW"]

println("Running algorithm $(alg)")

output = NDFW.traffic_assignment_lp(econfig)



if alg == "IFW" || alg == "CONIC" || alg == "BDM"
    @info "Running integer FW algorithm for Network Design"
    pv = missing
    penalty = missing
end



output_row = [
    ds,
    tn.number_of_nodes,
    tn.number_of_links,
    filealg,
    round(output.result, digits=sig_dig),
    output.base_cost,
    penalty,
    pv,
    fraction_removed,
    output.new_edges_opened,
    output.new_edge_count,
    output.cons_violation,
    round(output.total_time, digits=sig_dig),
    output.number_of_nodes,
    output.total_time / output.number_of_nodes,
    output.primal_dual_prog,
    output.abs_dual_gap,
    total_scenarios,
    solver
]

@info "Obtained solution:" round(output.result, digits=sig_dig)
push!(results, output_row)

CSV.write(results_file, results)

# CSV.write(
#     "$(folder_name)/solution_$(ds)_$(penalty)_$(pv)_$(filealg)_$(fraction_removed)_$(total_scenarios).csv", 
#     Tables.table(output.soln_vector),
#     writeheader=false
# )


