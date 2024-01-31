function build_conic_model(
    tn,
    r::Vector{Float64}, 
    num_nodes::Int, 
    edges::Vector{Tuple{Int,Int}}, 
    Z_set, 
    O_set, 
    demands, 
    removable_edges::Vector{Tuple{Int,Int}}, 
    max_flow::Vector{Float64},
    optimizer, 
    gparams,
    res
    )

    all_nodes = vcat(Z_set, O_set)
    S_set = collect(1:gparams.total_scenarios)
    #edges_with_cones = []

    m = Model()
    set_optimizer(m, optimizer)
    @variable(m, x[e in edges, z in Z_set, s in S_set] >= 0)
    @variable(m, y[i in removable_edges], Bin)
    #@variable(m, 0 <= y[i in removable_edges] <= 1)
    @variable(m, w[e in edges, s in S_set] >= 0)
    @variable(m, x_agg[e in edges, s in S_set])
    # w >= (∑_z x_ijz)ᵅ

    # flow balance constraints
    # O_set: nodes except for the zones
    # Z_set: zones
    @constraint(
        m,
        demand_balance[i in O_set, z in Z_set, s in S_set; i ∉ Z_set], 
        sum(x[(i, j), z, s] for j in all_nodes if (i, j) in edges) ==
        sum(x[(j, i), z, s] for j in all_nodes if (j, i) in edges)
    )

    #outflow constraints
    @constraint(
        m,
        outflow[iz in Z_set, oz in Z_set, s in S_set],
        sum(x[(iz, j), oz, s] for j in O_set if (iz, j) in edges) ==
        demands[s][oz][iz]
        )

    #inflow constraints
    @constraint(
        m,
        inflow[oz in Z_set, s in S_set],
        sum(x[(j, oz), oz, s] for j in O_set if (j, oz) in edges) ==
        sum(demands[s][oz][iz] for iz in Z_set)
        )

    # other network flow constraints
    @constraint(m, design_constraint[e in removable_edges, z in Z_set, s in S_set], x[e,z,s] <= max_flow[z] * y[e])
    

    obj_coefficients = get_obj_coefficients_cone(tn)
    obj_coefficients_forall_scenarios = repeat(obj_coefficients, gparams.total_scenarios)
    obj_coefficients_linear = get_obj_coefficients_linear(tn)
    obj_coefficients_const = get_obj_coefficients_constant(tn)
    obj_powers = get_obj_powers_cone(tn)

    @constraint(m,
    aggregate_destination_flows[e in edges, s in S_set], 
    x_agg[e, s] == sum(x[e, z, s] for z in Z_set)
    )
    
    

    @constraint(m,
        power_cone_cons_edge[e in edges, s in S_set],
        vcat(w[e, s], 1.0, x_agg[e, s]) in MOI.PowerCone(obj_powers[gparams.link_dic[e[1],e[2]]]),
    )

    # @constraint(m,
    #     power_cone_cons_edge[e in edges, s in S_set],
    #     vcat(w[e, s], 1.0, x_agg[e, s]) in MOI.Nonnegatives(3),
    # )
    
    # @constraint(m,
    #     power_cone_cons_edge[e in edges, s in S_set],
    #     vcat(w[e, s], 1.0, x_agg[e, s]) in MOI.SecondOrderCone(3),
    # )


    scale_factor = 1.0
    prob = 1.0/length(S_set)

    @objective(m, Min, 
        scale_factor * dot(r, y) +
        scale_factor * prob * sum(obj_coefficients[i] * w[edges[i], s] for i in 1:length(edges) for s in S_set) +
        scale_factor * prob * sum(obj_coefficients_linear[gparams.link_dic[e[1], e[2]]] * x_agg[e, s] for e in edges for s in S_set) +
        scale_factor * prob * sum(obj_coefficients_const)
        )

    return m, y, x, w
end

function get_obj_coefficients_cone(tn)
    obj_coeffcients = zeros(length(tn.init_node))

    for i in eachindex(obj_coeffcients)
        obj_coeffcients[i] += tn.free_flow_time[i] * (tn.b[i] / (tn.capacity[i]^tn.power[i]) / (tn.power[i] + 1))
        #obj_coeffcients[i] += tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
        #sum += 0.0001*x_agg[i]
    end

    return obj_coeffcients
end

function get_obj_coefficients_linear(tn)
    obj_coeffcients = zeros(length(tn.init_node))

    for i in eachindex(obj_coeffcients)
        obj_coeffcients[i] += tn.free_flow_time[i] 
    end

    return obj_coeffcients
end

function get_obj_coefficients_constant(tn)
    obj_coeffcients = zeros(length(tn.init_node))

    for i in eachindex(obj_coeffcients)
        obj_coeffcients[i] += tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
    end

    return obj_coeffcients
end

function get_obj_powers_cone(tn)
    obj_powers = zeros(length(tn.init_node))

    for i in eachindex(obj_powers)
        obj_powers[i] += 1/(tn.power[i] + 1)
    end

    return obj_powers
end

# function build_conic_perspective_model(
#     tn,
#     r::Vector{Float64},
#     num_nodes::Int,
#     edges::Vector{Tuple{Int,Int}},
#     Z_set,
#     O_set,
#     demands,
#     removable_edges::Vector{Tuple{Int,Int}},
#     max_flow::Vector{Float64},
#     optimizer,
#     alpha,
#     gparams
# )

#     all_nodes = vcat(Z_set, O_set)

#     m = Model()
#     set_optimizer(m, optimizer)
#     @variable(m, x[e in edges, z in Z_set] >= 0)
#     @variable(m, y[i in removable_edges], Bin)
#     @variable(m, τ[e in edges] >= 0)

#     #@variable(m, 0 <= y[i in removable_edges] <= 1)

#     # flow balance constraints
#     # O_set: nodes except for the zones
#     # Z_set: zones
#     @constraint(
#         m,
#         demand_balance[i in O_set, z in Z_set; i ∉ Z_set],
#         sum(x[(i, j), z] for j in all_nodes if (i, j) in edges) ==
#         sum(x[(j, i), z] for j in all_nodes if (j, i) in edges)
#     )

#     #outflow constraints
#     @constraint(
#         m,
#         outflow[iz in Z_set, oz in Z_set],
#         sum(x[(iz, j), oz] for j in O_set if (iz, j) in edges) ==
#         demands[oz][iz]
#     )

#     #inflow constraints
#     @constraint(
#         m,
#         inflow[oz in Z_set],
#         sum(x[(j, oz), oz] for j in O_set if (j, oz) in edges) ==
#         sum(demands[oz][iz] for iz in Z_set)
#     )

#     # other network flow constraints
#     #@constraint(m, design_constraint[e in removable_edges, z in Z_set], x[e, z] <= max_flow[z] * y[e])
#     # w >= (∑_z x_ijz)ᵅ

#     obj_coefficients = get_obj_coefficients_cone(tn)
#     obj_coefficients_linear = get_obj_coefficients_linear(tn)
#     obj_coefficients_const = get_obj_coefficients_constant(tn)
#     obj_powers = get_obj_powers_cone(tn)

#     @constraint(m,
#         power_cone_cons_remedge[e in removable_edges],
#         [τ[e], y[e], sum(x[e, z] for z in Z_set)] ∈ MOI.PowerCone(obj_powers[gparams.link_dic[e[1], e[2]]]),
#     )

#     println("adding power cone constraints with coefficient $(alpha)")
#     @constraint(m,
#         power_cone_cons_edge[e in edges; e ∉ removable_edges],
#         [τ[e], 1.0, sum(x[e, z] for z in Z_set)] ∈ MOI.PowerCone(obj_powers[gparams.link_dic[e[1], e[2]]]),
#     )

#     @objective(m, Min,
#         dot(r, y) +
#         dot(obj_coefficients, τ) +
#         sum(obj_coefficients_linear[gparams.link_dic[e[1], e[2]]] * sum(x[e, z] for z in Z_set) for e in edges) +
#         sum(obj_coefficients_const)
#     )

#     return m, y, x, τ
# end


function build_conic_perspective_model(
    tn,
    r::Vector{Float64}, 
    num_nodes::Int, 
    edges::Vector{Tuple{Int,Int}}, 
    Z_set, 
    O_set, 
    demands, 
    removable_edges::Vector{Tuple{Int,Int}}, 
    max_flow::Vector{Float64},
    optimizer, 
    gparams,
    res
    )

    all_nodes = vcat(Z_set, O_set)
    S_set = collect(1:gparams.total_scenarios)
    #edges_with_cones = edges[1:5]

    m = Model()
    set_optimizer(m, optimizer)
    @variable(m, x[e in edges, z in Z_set, s in S_set] >= 0)
    @variable(m, y[i in removable_edges], Bin)
    #@variable(m, 0 <= y[i in removable_edges] <= 1)
    @variable(m, w[e in edges, s in S_set] >= 0)
    @variable(m, x_agg[e in edges, s in S_set])
    # w >= (∑_z x_ijz)ᵅ

    # flow balance constraints
    # O_set: nodes except for the zones
    # Z_set: zones
    @constraint(
        m,
        demand_balance[i in O_set, z in Z_set, s in S_set; i ∉ Z_set], 
        sum(x[(i, j), z, s] for j in all_nodes if (i, j) in edges) ==
        sum(x[(j, i), z, s] for j in all_nodes if (j, i) in edges)
    )

    #outflow constraints
    @constraint(
        m,
        outflow[iz in Z_set, oz in Z_set, s in S_set],
        sum(x[(iz, j), oz, s] for j in O_set if (iz, j) in edges) ==
        demands[s][oz][iz]
        )

    #inflow constraints
    @constraint(
        m,
        inflow[oz in Z_set, s in S_set],
        sum(x[(j, oz), oz, s] for j in O_set if (j, oz) in edges) ==
        sum(demands[s][oz][iz] for iz in Z_set)
        )

    # other network flow constraints
    # @constraint(m, design_constraint[e in removable_edges, z in Z_set, s in S_set], x[e,z,s] <= max_flow[z] * y[e])
    

    obj_coefficients = get_obj_coefficients_cone(tn)
    obj_coefficients_forall_scenarios = repeat(obj_coefficients, gparams.total_scenarios)
    obj_coefficients_linear = get_obj_coefficients_linear(tn)
    obj_coefficients_const = get_obj_coefficients_constant(tn)
    obj_powers = get_obj_powers_cone(tn)

    @constraint(m,
    aggregate_destination_flows[e in edges, s in S_set], 
    x_agg[e, s] == sum(x[e, z, s] for z in Z_set)
    )
    
    @constraint(m,
        power_cone_cons_remedge[e in removable_edges, s in S_set],
        [w[e, s], y[e], x_agg[e, s]] ∈ MOI.PowerCone(obj_powers[gparams.link_dic[e[1], e[2]]]),
    )

    @constraint(m,
        power_cone_cons_edge[e in edges, s in S_set; e ∉ removable_edges],
        vcat(w[e, s], 1.0, x_agg[e, s]) in MOI.PowerCone(obj_powers[gparams.link_dic[e[1],e[2]]]),
    )

    # @constraint(m,
    #     power_cone_cons_edge_lin[e in edges, s in S_set; e ∉ edges_with_cones],
    #     w[e, s] >= x_agg[e, s]
    #     )

    scale_factor = 1.0
    prob = 1.0/length(S_set)

    @objective(m, Min, 
        scale_factor * dot(r, y) +
        scale_factor * prob * sum(obj_coefficients[i] * w[edges[i], s] for i in 1:length(edges) for s in S_set) +
        scale_factor * prob * sum(obj_coefficients_linear[gparams.link_dic[e[1], e[2]]] * x_agg[e, s] for e in edges for s in S_set) +
        scale_factor * prob * sum(obj_coefficients_const)
        )

    return m, y, x, w
end