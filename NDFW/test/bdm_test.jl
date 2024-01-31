using Test
using Boscia
using NDFW
using ForwardDiff
using LinearAlgebra
using Random
using StatsBase
using JuMP
using FrankWolfe
using ProfileView
``````
alg = "BDM"
#ds = "Berlin-Friedrichshain"
ds = "Braess-Example"
base_cost = 1000.0
rem_cd = true
total_scenarios = 5
scale = 0.1
fraction_removed = 0.03
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

optimizer, con_list = NDFW.initialize_model(alg, tn, g1, params, gparams, type_feas_cuts, total_scenarios)

lmo = NDFW.initialize_lmo(alg, optimizer, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)


optimizer2, con_list = NDFW.initialize_model("IFW", tn, g1, params, gparams, type_feas_cuts, total_scenarios)

o = optimizer2

NDFW.add_network_design_constraints!(
    o,
    new_edge_list,
    con_list,
    g1,
    params,
    total_scenarios
)

ilmo = NDFW.initialize_lmo("IFW", optimizer2, g1.graph, gparams, use_bigm, use_adaptive, use_reverse)



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

#x_nd = rand([0.0, 1.0], new_edge_count)
x_nd = zeros(new_edge_count)
x[var_count-new_edge_count+1:var_count] = x_nd


fbdm, gbdm = NDFW.build_f_and_grad("IFW", tn, gparams, cost_of_expansion, res, total_scenarios)

#auto_grad = ForwardDiff.gradient(f, x)


#grad = grad!(storage, x)
fb = fbdm(x)

storage = zeros(var_count)
gbdm(storage, x)
gb = copy(storage)


frac_to_bound = 0.2
vars_to_bound = Int(ceil(frac_to_bound * length(lmo.int_vars)))
selected_vars = StatsBase.sample(collect(1:length(lmo.int_vars)), vars_to_bound, replace=false)

for i in eachindex(lmo.int_vars)
    if i in selected_vars
        if x_nd[i] == 1.0
            sense = :greaterthan
        elseif x_nd[i] == 0.0
            sense = :lessthan
        else
            throw(ArgumentError("invalid value"))
        end

        Boscia.set_bound!(lmo, lmo.int_vars[i], x_nd[i], sense)
    end        
end

lmo_int_vars = Boscia.get_binary_variables(ilmo)
lower_bound_list = Boscia.get_lower_bound_list(ilmo)
upper_bound_list = Boscia.get_upper_bound_list(ilmo)



for i in eachindex(lmo_int_vars)
    if i in selected_vars
        rel_cons = 0
        if x_nd[i] == 1.0
            sense = :greaterthan

            for cx in lower_bound_list
                valid = Boscia.is_constraint_on_int_var(ilmo, cx, [lmo_int_vars[i].value])
                if valid == true
                    rel_cons = cx
                    break
                else
                    continue
                end
            end

        elseif x_nd[i] == 0.0
            sense = :lessthan

            for cx in upper_bound_list
                valid = Boscia.is_constraint_on_int_var(ilmo, cx, [lmo_int_vars[i].value])
                if valid == true
                    rel_cons = cx
                    break
                else
                    continue
                end
            end


        else
            throw(ArgumentError("invalid value"))
        end

        Boscia.set_bound!(ilmo, rel_cons, x_nd[i], sense)
    end
end



xb = Boscia.compute_extreme_point(lmo, gb)
xi = Boscia.compute_extreme_point(ilmo, gb)

println("objective for bdm $(fbdm(xb) - res)")
println("objective for ifw $(fbdm(xi) - res)")