using Test
using Boscia
using NDFW
using ForwardDiff
using LinearAlgebra
using Random
using StatsBase

alg = "BNDLMO-P"
ds = "Berlin-Friedrichshain"
#ds = "Braess-Example"
#base_cost = 1181.72068864639
base_cost = 77.2
rem_cd = true
total_scenarios = 4
scale = 0.1
fraction_removed = 0.1
seed = Int(round(1000 * (4 + fraction_removed)))
#res = 100
res = 618039.92/100
penalty = 1000.0
pv = 1.5
type_feas_cuts = ["point"]
use_bigm = false
use_adaptive = true

@info "Choosing random seed $(seed)"
#Random.seed!(seed)

tn, g1, params, gparams = NDFW.initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)

cost_of_expansion = []

@info "Choosing random seed $(seed)"
#Random.seed!(seed)

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


# o = optimizer

# NDFW.add_network_design_constraints!(
#     o,
#     new_edge_list,
#     con_list,
#     g1,
#     params,
#     total_scenarios
# )


# lmo = NDFW.initialize_lmo("IFW", optimizer, g1.graph, gparams, use_bigm, use_adaptive)

blmo = NDFW.initialize_lmo("BNDLMO-P", optimizer, g1.graph, gparams, use_bigm, use_adaptive)



var_count = total_scenarios * (od_pair_count * edge_count + edge_count) + new_edge_count
total_scenario_variables = od_pair_count * edge_count + edge_count

x = rand(var_count)

for scen in 1:total_scenarios
    first_var = (scen - 1) * total_scenario_variables

    for i in 1:edge_count
        curr_x_flows = [x[first_var+(z-1)*edge_count+i] for z in 1:od_pair_count]

        x[first_var+od_pair_count*edge_count+i] = sum(curr_x_flows)
    end
end

x_nd = rand([0.0, 1.0], new_edge_count)
#x_nd = ones(new_edge_count)
x[var_count-new_edge_count+1:var_count] = x_nd


# alg = "IFW"
# fifw, gifw =  NDFW.build_f_and_grad(alg, tn, gparams, cost_of_expansion, res, total_scenarios)

alg = "BNDLMO-P"
fblmo, gblmo = NDFW.build_f_and_grad_penalty_method(alg, tn, gparams, cost_of_expansion, res, total_scenarios, use_bigm, penalty, pv)

#auto_grad = ForwardDiff.gradient(f, x)


#grad = grad!(storage, x)
# fi = fifw(x)
fb = fblmo(x)

# storage = zeros(var_count)
# gifw(storage, x)
# gi = copy(storage)

storage = zeros(var_count)
gblmo(storage, x)
gb = copy(storage)

# lmo_int_vars = Boscia.get_binary_variables(lmo)
# lower_bound_list = Boscia.get_lower_bound_list(lmo)
# upper_bound_list = Boscia.get_upper_bound_list(lmo)


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
#     Boscia.set_bound!(lmo, rel_cons, x_nd[i], sense)
# end


# Boscia.free_model(lmo.o)

Boscia.check_feasibility(blmo)
xb = Boscia.compute_extreme_point(blmo, gb)
# xi = Boscia.compute_extreme_point(lmo, gi)

# @assert all(x_nd .== xi[end-new_edge_count+1:end])
@assert all(x_nd .== xb[end-new_edge_count+1:end])