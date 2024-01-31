using Test
using Boscia
using NDFW
using ForwardDiff
using LinearAlgebra
using Random
using StatsBase

alg = "BNDLMO-P"
ds = "Berlin-Tiergarten"
base_cost = 892.0247
rem_cd = true
total_scenarios = 1
scale = 0.1
fraction_removed = 0.01
seed = Int(round(1000 * (4 + fraction_removed)))
res = 100
penalty = 1000.0
pv = 1.5
use_bigm = true


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
x[var_count-new_edge_count+1:var_count] = rand([0, 1], new_edge_count)


f, grad! = NDFW.build_f_and_grad_penalty_method(alg, tn, gparams, cost_of_expansion, res, total_scenarios, use_bigm, penalty, pv)

auto_grad = ForwardDiff.gradient(f, x)

storage = zeros(var_count)

grad = grad!(storage, x)

println("gap between actual grad and autograd $(LinearAlgebra.norm(grad - auto_grad))")
