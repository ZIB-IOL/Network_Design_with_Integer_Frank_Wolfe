# Multicommodity graph with weights
mutable struct MC_graph_with_weights
    no_commodities::Int64
    graph::Graphs.SimpleDiGraph{Int64}
    node_demands_dict::Dict{Int64,Dict{Int64,Vector{Float64}}} #scenario -> dest -> origin
    edge_costs_dict::Dict{Tuple{Int64, Int64}, Vector{Float64}}
    edge_capacities_dict::Dict{Tuple{Int64, Int64}, Vector{Float64}}
    edge_cum_cap_dict::Dict{Tuple{Int64, Int64}, Float64}
    number_of_zones::Int64
end

# Multicommodity graph with weights
mutable struct Params
    no_edges::Int64
    no_nodes::Int64
    no_zones::Int64
    edge_list
    edge_dict
    incoming_edges
    outgoing_edges
    node_demands_dict
    MAX_FLOW
end

mutable struct GraphParams
    init_nodes::Vector{Int64}
    term_nodes::Vector{Int64}
    travel_demand::Dict{Int64, Matrix{Float64}}
    first_thru_node::Int64
    link_dic
    od_pair_count::Int64
    removed_edges
    removed_edges_costs::Vector{Float64}
    MAX_FLOW::Float64
    total_scenarios::Int64
end

#Initialize empty constructor
GraphParams() = GraphParams([], [], Dict{Int64,Matrix{Float64}}(), 0, [], 0, [], [], Inf, 1)

struct Solution
    ta_data_name::String
    alg::String
    result::Float64
    base_cost::Float64
    new_edges_opened::Int
    new_edge_count::Int
    cons_violation::Float64
    soln_vector::Vector{Float64}
    new_edge_list::Vector{Tuple{Int64, Int64}}
    number_of_nodes::Int64
    total_time::Float64
    primal_dual_prog::Float64
    abs_dual_gap::Float64
end

Solution(
    ta_data_name::String, 
    alg::String,
    result::Float64,
    base_cost::Float64,
    new_edges_opened::Int,
    new_edge_count::Int,
    soln_vector::Vector{Float64}) = Solution(
        ta_data_name,
        alg,
        result,
        base_cost,
        new_edges_opened,
        new_edge_count,
        Inf,
        soln_vector,
        [],
        0,
        0.0,
        Inf,
        Inf)

struct ExperimentConfig
    ta_data_name::String
    alg::String
    normalize::Bool
    time_limit::Int64
    penalty::Float64
    max_fw_iter::Int64
    pv::Float64
    sig_dig::Int64
    remove_circular_demand::Bool
    type_feas_cuts::Vector{String}
    fraction_removed::Float64
    seed::Int64
    rel_dual_gap::Float64
    bs_verb::Bool
    fw_verb::Bool
    ls_verb::Bool
    dual_tightening::Bool
    total_scenarios::Int64
    use_bigm::Bool
    use_adaptive::Bool
    use_reverse::Bool
    scale::Float64
    type_of_nd_constraints::String
    solver::String
end

