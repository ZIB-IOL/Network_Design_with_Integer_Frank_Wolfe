using Test
using Boscia
using NDFW
using ForwardDiff
using LinearAlgebra
using Random
using StatsBase
using CSV
using DataFrames

alg = "IFW"
ds = "Berlin-Friedrichshain"
#ds = "Braess-Example"
base_cost = 1181.72068864639
#base_cost = 77.2
rem_cd = true
total_scenarios = 1
scale = 0.1
fraction_removed = 0.01
seed = Int(round(1000 * (4 + fraction_removed)))
#res = 100
res = 618039.92/100
penalty = 1000.0
pv = 1.5
type_feas_cuts = ["point"]

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


lmo = NDFW.initialize_lmo("IFW", optimizer, g1.graph, gparams)

blmo = NDFW.initialize_lmo("BNDLMO-P", optimizer, g1.graph, gparams)



var_count = total_scenarios * (od_pair_count * edge_count + edge_count) + new_edge_count

println(pwd())
yi = CSV.read("./examples/solution_Berlin-Friedrichshain_missing_missing_IFW_0.01_1.csv", DataFrame, header=false)[:, 1]
yb = CSV.read("./examples/blmo_test_v4/solution_Berlin-Friedrichshain_1000.0_1.5_BNDLMO-P_0.01_1.csv", DataFrame, header=false)[:, 1]

alg = "IFW"
fifw, gifw =  NDFW.build_f_and_grad(alg, tn, gparams, cost_of_expansion, res, total_scenarios)

alg = "BNDLMO-P"
fblmo, gblmo = NDFW.build_f_and_grad_penalty_method(alg, tn, gparams, cost_of_expansion, res, total_scenarios, penalty, pv)

totsum = 0.0

removed_edge_idx = [gparams.link_dic[removed_edge[1], removed_edge[2]] for removed_edge in gparams.removed_edges]

td_2_dest = [maximum([sum(gparams.travel_demand[s][:, z]) for s in 1:total_scenarios]) for z in 1:od_pair_count]

x_nd = yb[end - new_edge_count + 1:end]

total_scenario_vars = od_pair_count * edge_count + edge_count

# totsum = 0.0

# for scen in 1:total_scenarios
#     first_var = (scen - 1) * total_scenario_vars
#     #curr_x = @view(x[(first_var+1):(first_var+od_pair_count*edge_count)])

#     for j in 1:new_edge_count
#         for z in 1:od_pair_count
#             curr_x = @view(yb[(first_var+(z-1)*edge_count+1):(first_var+z*edge_count)])
#             #x_agg = 
#             #dest_edge_loc = (z - 1) * edge_count + removed_edge_idx[j]
#             totsum += penalty * max(curr_x[removed_edge_idx[j]] - td_2_dest[z] * x_nd[j], 0)^pv
#         end
#     end
# end

# bigM = tn.total_od_flow

# violation = []
# for scen = 1:total_scenarios
#     first_var = (scen - 1) * total_scenario_vars
#     x_agg = yb[(first_var+od_pair_count*edge_count+1):(first_var+od_pair_count*edge_count+edge_count)]
#     cons_violation = [max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0) for j in eachindex(x_nd)]
#     append!(violation, penalty * sum(cons_violation))
# end

# avg_violation = sum(violation) / length(violation)


#auto_grad = ForwardDiff.gradient(f, x)

#Check actual constraint violation


#grad = grad!(storage, x)
fi = fifw(yi)
fb = fblmo(yb)

storage = zeros(var_count)
gifw(storage, yi)
gi = copy(storage)

storage = zeros(var_count)
gblmo(storage, yb)
global gb = copy(storage)

lmo_int_vars = Boscia.get_binary_variables(lmo)
lower_bound_list = Boscia.get_lower_bound_list(lmo)
upper_bound_list = Boscia.get_upper_bound_list(lmo)


for i in eachindex(blmo.int_vars)
    if x_nd[i] == 1.0
        sense = :greaterthan
    elseif x_nd[i] == 0.0
        sense = :lessthan
    else
        throw(ArgumentError("invalid value"))
    end

    Boscia.set_bound!(blmo, blmo.int_vars[i], x_nd[i], sense)
end

Boscia.set_bound!(blmo, blmo.int_vars[4], 0.0, :greaterthan)

# for i in eachindex(lmo_int_vars)
#     rel_cons = 0
#     if x_nd[i] == 1.0
#         sense = :greaterthan

#         for cx in lower_bound_list
#             valid = Boscia.is_constraint_on_int_var(lmo, cx, [lmo_int_vars[i].value])
#             if valid == true
#                 rel_cons  = cx
#                 break
#             else
#                 continue
#             end
#         end

#     elseif x_nd[i] == 0.0
#         sense = :lessthan

#         for cx in upper_bound_list
#             valid = Boscia.is_constraint_on_int_var(lmo, cx, [lmo_int_vars[i].value])
#             if valid == true
#                 rel_cons = cx
#                 break
#             else
#                 continue
#             end
#         end


#     else
#         throw(ArgumentError("invalid value"))
#     end

#     println(rel_cons)
#     #Boscia.set_bound!(lmo, rel_cons, x_nd[i], sense)
#     Boscia.set_bound!(lmo, rel_cons, 1.0, :greaterthan)
# end


for i in eachindex(lmo_int_vars)
    rel_cons = 0

    for cx in lower_bound_list
        valid = Boscia.is_constraint_on_int_var(lmo, cx, [lmo_int_vars[i].value])
        if valid == true
            rel_cons = cx
            break
        else
            continue
        end
    end

    println(rel_cons)
    Boscia.set_bound!(lmo, rel_cons, 1.0, :greaterthan)
end

Boscia.free_model(lmo.o)

vbb = Boscia.compute_extreme_point(blmo, gb)
vbi = Boscia.compute_extreme_point(lmo, gb)

println("difference in solutions: $(norm(vbb[1:12552] - vbi[1:12552]))")

diff = vbb[1:12552] - vbi[1:12552]

diffnz = abs.(diff) .> 1e-10

varnz = collect(1:12552)[diffnz]

rel_edges = [e%523 for e in varnz]

rel_zones = [ceil(e/523) for e in varnz]

# @assert all(x_nd .== xi[end-new_edge_count+1:end])
# @assert all(x_nd .== xb[end-new_edge_count+1:end])