using Test
using Boscia
using NDFW
using ForwardDiff
using LinearAlgebra
using Random
using StatsBase
using JuMP
using FrankWolfe
using Gurobi
using Graphs
``````
alg = "BDM"
ds = "Berlin-Friedrichshain"
#ds = "Braess-Example"
base_cost = 1000.0
rem_cd = true
total_scenarios = 1
scale = 0.1
fraction_removed = 0.1
seed = Int(round(1000 * (4 + fraction_removed))) + 4
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
rel_dual_gap_ifw = 1e-3

@info "Choosing random seed $(seed)"
Random.seed!(seed)

tn, g1, params, gparams = NDFW.initialize_parameters(ds, alg, base_cost, rem_cd, total_scenarios, scale)


cost_of_expansion = []

@info "Choosing random seed $(seed)"
Random.seed!(seed)

edge_count = length(gparams.init_nodes)
cost_vector = 1.0 .+ rand(edge_count)


#arcs_to_be_expanded = [rand() > 0.5 for i in 1:length(edge_list)]
count_arcs_to_remove = Int(ceil(fraction_removed * length(params.edge_list)))
new_edge_list = StatsBase.sample(params.edge_list, count_arcs_to_remove, replace=false)
new_edge_list = [(47,48)]
new_edge_count = length(new_edge_list)


@info "Removed $(new_edge_count) our of $(length(params.edge_list)) edges"
cost_of_expansion = base_cost * ones(new_edge_count)

gparams.removed_edges = new_edge_list
gparams.removed_edges_costs = cost_of_expansion

od_pair_count = gparams.od_pair_count


var_count = total_scenarios * (od_pair_count * edge_count + edge_count) + new_edge_count
total_scenario_variables = od_pair_count * edge_count + edge_count

edge_list = [(tn.init_node[i], tn.term_node[i]) for i in 1:length(tn.init_node)]
y_val = rand([0.0, 1.0],new_edge_count)

GRB_ENV = Gurobi.Env()
travel_time_all = [zeros(od_pair_count * edge_count);cost_vector]

nc = Graphs.nv(g1.graph)
ie = length(edge_list)
distmx = Inf * ones(nc, nc)

scen = 1

st, dual_var_r, dual_var_t, dual_var_s, q = NDFW.all_or_nothing_LP(scen, travel_time_all, gparams, y_val, nc, edge_list, GRB_ENV)

temp_graph = copy(g1.graph)
modified_graph, curr_removed_edges = NDFW.modify_graph(temp_graph, y_val, gparams.removed_edges)
modified_costs = NDFW.modify_cost_vector(travel_time_all, gparams, curr_removed_edges, ie)
new_edge_list = [(edge.src, edge.dst) for edge in edges(modified_graph)]

orig_v, r, t, rp, p = NDFW.all_or_nothing_stochastic_detailed(modified_costs, gparams, modified_graph, g1.graph, scen, distmx)

v = NDFW.map_to_original(orig_v, temp_graph, gparams, total_scenarios)

println("shortest path solution_agg: $(dot(v[od_pair_count * edge_count + 1:end], cost_vector))")
println("shortest path solution_sep: $(dot(v[1:od_pair_count * edge_count], travel_time_all[1:od_pair_count * edge_count]))")


all_or_nothing_costs = 0.0

prob = 1/total_scenarios

for scen in 1:total_scenarios
    global all_or_nothing_costs += prob * sum(
        sum((dual_var_r[i, z] - dual_var_t[z]) * gparams.travel_demand[scen][i, z] for i in 1:od_pair_count)
        for z in 1:od_pair_count
    )
end

println("demand costs: $(all_or_nothing_costs)")