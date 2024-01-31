"""
The goal of this script is to compare the performance of the different algorithms specifically the biconjugate FW from the TA Librart and 
the BPCG methods (along with others) from the FrankWolfe.jl library. 
"""

using DataFrames
using CSV
using NDFW
using Dates
using Tables
using Plots


# ds = ARGS[1]
# alg = ARGS[2]
# solver = ARGS[3]

# println(ARGS[1])
# println(ARGS[2])
# println(ARGS[3])




ds = "Berlin-Friedrichshain"
alg = "SP"
solver="SCIP"


use_adaptive = true
type_of_nd_constraints = ""
filealg = alg
total_scenarios = 1
use_bigm = false
use_reverse = false



println("Running dataset $(ds)")

time_limit = 10800
normalize = true
#penalties = [1e1, 1e2, 1e3]
max_fw_iter = 2000
#powers = [1.5]

sig_dig = 3
exp_id = 6
rem_cd = true
type_feas_cuts = ["point"]
seed = 1
rel_dual_gap = 1e-3
bs_verb = true
fw_verb = true
ls_verb = false
dual_tightening = true
#total_scenarios = 12
scale = 0.1
norm_const = 1.0


econfig = NDFW.ExperimentConfig(
    ds,
    alg,
    normalize,
    time_limit,
    0.0,
    max_fw_iter,
    0.0,
    sig_dig,
    rem_cd,
    type_feas_cuts,
    0.0,
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

folder_name = "ta_experiment_afw"
#folder_name = "experiment_$(Dates.format(date_today, "ddmmyy"))_$(exp_id)" 

isdir(folder_name) || mkdir(folder_name)

results_file = "$(folder_name)/results_$(ds)_$(filealg)_$(total_scenarios)_$(solver).csv"

@info "Evaluating data set $(ds)"
tn = NDFW.load_ta_network(ds)

@info "Running tailored FW algorithm"
t = @elapsed _, _, res, ta_traj = NDFW.modified_ta_frank_wolfe(tn, log=:on, norm_const=norm_const)
#alg = "TN"
@info "Obtained solution:" round(res,digits=sig_dig)
output_row = [ds, 
                tn.number_of_nodes, 
                tn.number_of_links, 
                "TN", 
                round(res,digits=sig_dig), 
                missing, 
                missing, 
                missing, 
                missing, 
                missing, 
                missing, 
                missing,
                round(t,digits=sig_dig),
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

if normalize == true
    sres = res #used for scaling objective
    cost_of_expansion = round((res / tn.number_of_links), digits = 5) #cost of expansion set to average cost of flow across arcs
else
    sres = 100
end

println("Running algorithm $(alg)")

if alg == "LP"
    output = NDFW.traffic_assignment_lp(econfig)
elseif alg == "SP"
    t = @elapsed output, bpcg_traj = NDFW.traffic_assignment_sp(econfig, norm_const)
else
    throw(ArgumentError("Invalid alg argument $(alg)"))
end


output_row = [
    ds, 
    tn.number_of_nodes, 
    tn.number_of_links, 
    filealg, 
    round(output.result ,digits=sig_dig), 
    output.base_cost, 
    missing,
    missing,
    0,
    output.new_edges_opened, 
    output.new_edge_count, 
    output.cons_violation, 
    round(output.total_time, digits=sig_dig),
    output.number_of_nodes,
    output.total_time / output.number_of_nodes ,
    output.primal_dual_prog,
    output.abs_dual_gap,
    total_scenarios,
    solver
    ]

@info "Obtained solution:" round(output.result, digits=sig_dig)
push!(results, output_row)

#CSV.write(results_file, results)


ta_obj = [st.obj for st in ta_traj]
bpcg_obj = [st[2] for st in bpcg_traj]

plot(bpcg_obj, xaxis=:log10, yaxis=:log10, label="BPCG", linewidth=2)
plot!(ta_obj, xaxis=:log10, yaxis=:log10, label="BFW", linewidth=2)
title!("Comparison of FW methods")
xlabel!("Iterations")
ylabel!("Objective Value")

savefig("$(ds)_objective_comparison.pdf")


ta_gap = [st.dual_gap for st in ta_traj]
bpcg_gap = [st[4] for st in bpcg_traj]

plot(bpcg_gap, xaxis=:log10, yaxis=:log10, label="BPCG", linewidth=2)
plot!(ta_gap, xaxis=:log10, yaxis=:log10, label="BFW", linewidth=2)
title!("Comparison of FW methods")
xlabel!("Iterations")
ylabel!("FW Gap")

savefig("$(ds)_gap_comparison.pdf")

