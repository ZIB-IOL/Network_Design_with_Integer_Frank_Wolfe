#construct dictionary for incoming and outgoing edges for each vertex in graph
function neighbouring_edges(g)
    incoming_edges = Dict{Int64, Vector{Int64}}()
    outgoing_edges = Dict{Int64, Vector{Int64}}()

    for (i, edge) in enumerate(edges(g))
        src = edge.src
        dst = edge.dst

        if haskey(incoming_edges, dst)
            append!(incoming_edges[dst], i)
        else
            incoming_edges[dst] = [i]
        end

        if haskey(outgoing_edges, src)
            append!(outgoing_edges[src], i)
        else
            outgoing_edges[src] = [i]
        end
    end
    return incoming_edges, outgoing_edges
end

function add_mc_edge_wts!(g::MC_graph_with_weights, src, dst, no_commodities, weight, capacity, cum_capacity)
    Graphs.add_edge!(g.graph, src, dst)
    g.edge_costs_dict[(src, dst)] = zeros(no_commodities)
    g.edge_capacities_dict[(src, dst)] = zeros(no_commodities)
    for comm in range(1, step=1, stop=g.no_commodities)
        g.edge_costs_dict[(src, dst)][comm] = weight[comm]
        g.edge_capacities_dict[(src, dst)][comm] = capacity[comm]
    end
    if cum_capacity == NaN
        g.edge_cum_cap_dict[(src, dst)] = sum(g.edge_capacities_dict[(src, dst)])
    else
        g.edge_cum_cap_dict[(src, dst)] = cum_capacity
    end
end


function create_mc_cost_coeffs(g::MC_graph_with_weights, cost_of_expansion)
    edge_list = [(edge.src, edge.dst) for edge in collect(edges(g1.graph))]
    cost_coeffs = []
    for comm in range(1, step=1, stop=g.no_commodities)
        for edge in edge_list
            append!(cost_coeffs, g.edge_costs_dict[edge][comm])
        end
    end
    #append!(cost_coeffs, cost_of_expansion*ones(no_new_edges))

    return cost_coeffs
end

"""
returns the variable location given destination and edge
"""
function vl(dest, src, dst, params)
    edge = params.edge_dict[(src, dst)]
    return (dest - 1)*params.no_edges + edge
end

function revloc(loc, params)
    zone = Int64(ceil(loc/params.no_edges))
    edge = loc%params.no_edges == 0 ? params.no_edges : loc%params.no_edges
    return zone, params.edge_list[edge]
end

function build_constraint(type, curr_constraint, dest_zone, edge_list, iterator)
    if type == "inflow"
        for i in iterator
            curr_edge = edge_list[i]
            if curr_constraint == ""
                curr_constraint = curr_constraint*"-x_$(dest_zone),$(curr_edge[1])a$(curr_edge[2])"
            else
                curr_constraint = curr_constraint*" -x_$(dest_zone),$(curr_edge[1])a$(curr_edge[2])"
            end
        end
    elseif type == "outflow"
        for i in iterator
            curr_edge = edge_list[i]
            if curr_constraint == ""
                curr_constraint = curr_constraint*"x_$(dest_zone),$(curr_edge[1])a$(curr_edge[2])"
            else
                curr_constraint = curr_constraint*"+x_$(dest_zone),$(curr_edge[1])a$(curr_edge[2])"
            end
        end
    end
    return curr_constraint
end


function remove_circular_demand(tn)
    """
    Checkis if there is any circular demand in the network and removes it
    """
    for source in 1:size(tn.travel_demand)[1]
        for dest in 1:size(tn.travel_demand)[2]
            if (source == dest) & (tn.travel_demand[source, dest] > 0)
                @info "Removing circular demand from $(source) to $(dest)"
                tn.travel_demand[source, dest] = 0
            end
        end
    end

    return tn
end

