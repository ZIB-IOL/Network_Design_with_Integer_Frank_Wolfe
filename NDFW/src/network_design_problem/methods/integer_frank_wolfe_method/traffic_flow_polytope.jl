"""
Flow polytope without zones. This is to be used when the number of nodes are equal to number of zones. 
which indicates that there are no zones. 
INPUT:
optimizer: Optimizer with model
g1: graph data structures
params: set of graph parameters
OUTPUT:
optimizer: updated optimizer with new constraints
constraint_dict: set of constraints added
"""
function create_flow_polytope(
    optimizer,
    g1::MC_graph_with_weights,
    params::Params,
    scenario_cnt::Int64
    )

    noe = params.no_edges
    non = params.no_nodes
    var_count = params.no_edges * non
    incoming_edges = params.incoming_edges
    outgoing_edges = params.outgoing_edges
    edge_list = params.edge_list

    var_count = noe * non

    #make sure to silence hte optimizer output
    MOI.set(optimizer, MOI.Silent(), true)

    # I proceed with the assumption that the indices are commodity wise. 
    x = MOI.add_variables(optimizer, var_count)
    x_agg = MOI.add_variables(optimizer, noe)

    #y = MOI.add_variable(optimizer)
    #MOI.add_constraint(optimizer, y, MOI.ZeroOne())

    constraint_dict = Dict()
    constraint_dict["flow"] = []

    for dest in range(1, step=1, stop=non)
        for node in range(1, step=1, stop=non)
            #Initialize incoming edges with weight = 1
            if node != dest
                inc_terms = []
                out_terms = []

                curr_constraint = ""

                if haskey(incoming_edges, node)
                    inc_terms = [MOI.ScalarAffineTerm(-1.0, x[(dest - 1)*noe + i]) for i in incoming_edges[node]]

                    for i in incoming_edges[node]
                        curr_edge = edge_list[i]
                        if curr_constraint == ""
                            curr_constraint = curr_constraint*"- x_$(dest),$(curr_edge[1]) a $(curr_edge[2])"
                        else
                            curr_constraint = curr_constraint*" - x_$(dest),$(curr_edge[1]) a $(curr_edge[2])"
                        end
                    end
                end
                if haskey(outgoing_edges, node)
                    out_terms = [MOI.ScalarAffineTerm(1.0, x[(dest - 1)*noe + i]) for i in outgoing_edges[node]]
                    for i in outgoing_edges[node]
                        curr_edge = edge_list[i]
                        if curr_constraint == ""
                            curr_constraint = curr_constraint*"+ x_$(dest),$(curr_edge[1]) a $(curr_edge[2])"
                        else
                            curr_constraint = curr_constraint*" + x_$(dest),$(curr_edge[1]) a $(curr_edge[2])"
                        end
                    end
                end

                inc_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(inc_terms)
                out_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(out_terms)

                rhs = g1.node_demands_dict[scenario_cnt][dest][node]

                MOI.add_constraint(
                    optimizer, 
                    MOI.ScalarAffineFunction([inc_terms; out_terms], 0.0),
                    MOI.EqualTo(rhs)
                    )

                curr_constraint = curr_constraint*" == $(g1.node_demands_dict[dest][node])"
                append!(constraint_dict["flow"], [curr_constraint])
            end
        end
    end

    for edge in range(1, step=1, stop=noe)
        agg_term = [MOI.ScalarAffineTerm(1.0, x_agg[edge])]
        split_terms = [MOI.ScalarAffineTerm(-1.0, x[(dest - 1)*noe + edge]) for dest in range(1, step=1, length=non)]

        agg_term = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(agg_term)
        split_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(split_terms)

        MOI.add_constraint(
                optimizer, 
                MOI.ScalarAffineFunction([agg_term; split_terms], 0.0),
                MOI.EqualTo(0.0)
                )
    end

    for dest in range(1, step=1, stop=non)
        for edge in range(1, step=1, stop=noe)
            MOI.add_constraint(optimizer, x[(dest - 1)*noe + edge], MOI.GreaterThan(0.0))
            #MOI.add_constraint(optimizer, x[(dest - 1)*no_edges + edge], MOI.Integer())
        end
    end

    return optimizer, constraint_dict
end



"""
Flow polytope with zones. This is to be used when the number of nodes are different then the number of zones. 
which indicates that there are no zones. This function is primarily initialized the variables and then calls 
other functions which add each type of constraint.  
INPUT:
optimizer: Optimizer with model
g1: graph data structures
params: set of graph parameters
OUTPUT:
optimizer: updated optimizer with new constraints
constraint_dict: set of constraints added
"""
function create_flow_polytope_zones(
    optimizer,
    g1::MC_graph_with_weights, 
    params::Params,
    scenario_cnt::Int64
    )

    noe = params.no_edges
    noz = params.no_zones
    non = params.no_nodes
    edge_list = params.edge_list
    var_count = params.no_edges * noz
    incoming_edges = params.incoming_edges
    outgoing_edges = params.outgoing_edges
    node_demands_dict = params.node_demands_dict

    #make sure to silence hte optimizer output
    #MOI.set(optimizer, MOI.Silent(), true)

    x = MOI.add_variables(optimizer, var_count)  # Variables for flow on each arc to each destination
    x_agg = MOI.add_variables(optimizer, noe) #Aggregate flow on each arc

    #Name the variables. Useful when model is output. 
    for i in 1:var_count
        zone, (src, dst) = revloc(i, params)
        MOI.set(optimizer, MOI.VariableName(), x[i], "x_$(scenario_cnt),$(zone),$(src) a $(dst)")
    end

    #@info "Initialized $(var_count) variables from $(x[1]) to $(x[end])"

    for i in 1:noe
        zone, (src, dst) = revloc(i, params)
        MOI.set(optimizer, MOI.VariableName(), x_agg[i], "xagg_$(scenario_cnt),$(zone),$(src) a $(dst)")
    end

    #@info "Initialized $(noe) variables from $(x_agg[1]) to $(x_agg[end])"

    #Dictionary of constraints 
    constraint_dict = Dict()
    constraint_dict["flow"] = []
    constraint_dict["zoneoutflow"] = []
    constraint_dict["zoneoutflow_c"] = []
    constraint_dict["zoneinflow"] = []
    constraint_dict["aggregation"] = []

    constraint_dict["network_design"] = []

    #Add flow constraints for nodes not in the zones. 
    for dest_zone in range(1, step=1, stop=noz)
        for node in range(1, step=1, stop=non)
            #Initialize incoming edges with weight = 1
            if (!(node in 1:noz)) & (haskey(outgoing_edges, node) | haskey(incoming_edges, node))
                add_flow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
            end
        end
    end

    #@info "Initialized $(length(constraint_dict["flow"])) flow constraints for nodes except for zones"

    #Add outflow constraints for nodes in the zones
    for dest_zone in range(1, step=1, stop=noz)
        for node in range(1, step=1, stop=noz)
            if dest_zone != node
                if haskey(outgoing_edges, node)
                    add_outflow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
                end
            end
        end
    end

    #@info "Initialized $(length(constraint_dict["zoneoutflow"])) outflow constraints for zones with destinations different from themselves"

    #Add inflow constraints for nodes in the zones. 
    for node in range(1, step=1, stop=noz)
        dest_zone = node
        if haskey(incoming_edges, node)
            add_inflow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
        end
    end

    #@info "Initialized $(length(constraint_dict["zoneinflow"])) inflow constraints for zones"


    #Add constraints to manage or limit circular flows. 
    for node in range(1, step=1, stop=noz)
        dest_zone = node
        if haskey(outgoing_edges, node)
            add_zoneoutflow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
        end
    end

    #@info "Initialized $(length(constraint_dict["zoneoutflow_c"])) outflow constraints for zones with destinations back to themselves"

    #Add constraints to aggregate the flows to different destinatios on each arc
    for edge in range(1, step=1, stop=noe)
        add_aggregation_constraint(optimizer, x, x_agg, edge, constraint_dict, g1, params)
    end

    #@info "Initialized $(length(constraint_dict["aggregation"])) aggregation constraints for edges"

    #Add non nonnegative constraints for arc flows
    for idx in 1:var_count
        MOI.add_constraint(optimizer, x[idx], MOI.GreaterThan(0.0))
    end

    #@info "Initialized $(var_count) nonnegative constraints for edges"

    #Add non nonnegative constraints for aggregated arc flows
    for idx in 1:noe
        MOI.add_constraint(optimizer, x_agg[idx], MOI.GreaterThan(0.0))
    end

    #@info "Initialized $(noe) nonnegative constraints for aggregate edges"

    # add_network_design_constraints_internally!(
    # optimizer,
    # x_agg, 
    # y,
    # new_edge_list,
    # constraint_dict,
    # g1, 
    # params::Params)

    return optimizer, constraint_dict
end

"""
Add flow constraint which ensure that inflow = outflow for all nodes which are not zones. 
INPUT:
optimizer: optimization model
x : variable
dest_zone : destination of specified flow
node: node on which flow is being constrained
constraint_dict: dictionary of all constraints
g1: graph data strcuture
params: other graph parameters
OUTPUT:
None
"""
function add_flow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
    inc_terms = []
    out_terms = []

    curr_constraint = ""

    edge_list = params.edge_list
    incoming_edges = params.incoming_edges
    outgoing_edges = params.outgoing_edges
    node_demands_dict = params.node_demands_dict

    if haskey(incoming_edges, node)
        inc_terms = [MOI.ScalarAffineTerm(-1.0, x[vl(dest_zone, edge_list[i][1], edge_list[i][2], params)]) for i in incoming_edges[node]]

        curr_constraint = build_constraint("inflow", curr_constraint, dest_zone, edge_list, incoming_edges[node])
    end
    if haskey(outgoing_edges, node)
        out_terms = [MOI.ScalarAffineTerm(1.0, x[vl(dest_zone, edge_list[i][1], edge_list[i][2], params)]) for i in outgoing_edges[node]]

        curr_constraint = build_constraint("outflow", curr_constraint, dest_zone, edge_list, outgoing_edges[node])
    end

    inc_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(inc_terms)
    out_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(out_terms)

    rhs = 0.0

    c = MOI.add_constraint(
        optimizer, 
        MOI.ScalarAffineFunction([inc_terms; out_terms], 0.0),
        MOI.EqualTo(rhs)
        )

    curr_constraint = "f" * curr_constraint * "==$(g1.node_demands_dict[scenario_cnt][dest_zone][node])"
    MOI.set(optimizer, MOI.ConstraintName(), c, curr_constraint)

    append!(constraint_dict["flow"], [curr_constraint])
end

"""
Add outflow constraint which ensure that outflow = supply for all nodes which zones. 
INPUT:
optimizer: optimization model
x : variable
dest_zone : destination of specified flow
node: node on which flow is being constrained
constraint_dict: dictionary of all constraints
g1: graph data strcuture
params: other graph parameters
OUTPUT:
None
"""
function add_outflow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
    out_terms = []

    edge_list = params.edge_list
    outgoing_edges = params.outgoing_edges
    node_demands_dict = params.node_demands_dict

    curr_constraint = ""
    out_terms = [MOI.ScalarAffineTerm(1.0, x[vl(dest_zone, edge_list[i][1], edge_list[i][2], params)]) for i in outgoing_edges[node]]


    curr_constraint = build_constraint("outflow", curr_constraint, dest_zone, edge_list, outgoing_edges[node])

    out_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(out_terms)

    rhs = node_demands_dict[scenario_cnt][dest_zone][node]

    c = MOI.add_constraint(
        optimizer, 
        MOI.ScalarAffineFunction(out_terms, 0.0),
        MOI.EqualTo(rhs)
        )
    
    curr_constraint = "o" * curr_constraint * "==$(g1.node_demands_dict[scenario_cnt][dest_zone][node])"

    MOI.set(optimizer, MOI.ConstraintName(), c, curr_constraint)

    append!(constraint_dict["zoneoutflow"], [curr_constraint])
end

"""
Add flow constraint which ensure that inflow = - total supply at current zone, for all nodes which are zones. 
INPUT:
optimizer: optimization model
x : variable
dest_zone : destination of specified flow
node: node on which flow is being constrained
constraint_dict: dictionary of all constraints
g1: graph data strcuture
params: other graph parameters
OUTPUT:
None
"""
function add_inflow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
    inc_terms = []

    incoming_edges = params.incoming_edges
    edge_list = params.edge_list
    node_demands_dict = params.node_demands_dict
    curr_constraint = ""

    inc_terms = [MOI.ScalarAffineTerm(-1.0, x[vl(dest_zone, edge_list[i][1], edge_list[i][2], params)]) for i in incoming_edges[node]]


    curr_constraint = build_constraint("inflow", curr_constraint, dest_zone, edge_list, incoming_edges[node])

    inc_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(inc_terms)

    rhs = -sum(node_demands_dict[scenario_cnt][dest_zone])

    c = MOI.add_constraint(
        optimizer, 
        MOI.ScalarAffineFunction(inc_terms, 0.0),
        MOI.EqualTo(rhs)
        )

    curr_constraint = "i"*curr_constraint*"==$(rhs)"

    MOI.set(optimizer, MOI.ConstraintName(), c, curr_constraint)

    append!(constraint_dict["zoneinflow"], [curr_constraint])
end

"""
Add flow constraint for all flows which consider the origin as destination. 
INPUT:
optimizer: optimization model
x : variable
dest_zone : destination of specified flow
node: node on which flow is being constrained
constraint_dict: dictionary of all constraints
g1: graph data strcuture
params: other graph parameters
OUTPUT:
None
"""
function add_zoneoutflow_constraint!(optimizer, x, dest_zone, node, constraint_dict, g1, params, scenario_cnt)
    out_terms = []
    edge_list = params.edge_list
    outgoing_edges = params.outgoing_edges
    node_demands_dict = params.node_demands_dict

    curr_constraint = ""

    out_terms = [MOI.ScalarAffineTerm(1.0, x[vl(dest_zone, edge_list[i][1], edge_list[i][2], params)]) for i in outgoing_edges[node]]


    curr_constraint = build_constraint("outflow", curr_constraint, dest_zone, edge_list, outgoing_edges[node])

    out_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(out_terms)

    dem = node_demands_dict[scenario_cnt][dest_zone][node]
    rhs = dem

    c= MOI.add_constraint(
        optimizer, 
        MOI.ScalarAffineFunction(out_terms, 0.0),
        MOI.EqualTo(rhs)
        )

    curr_constraint = "zo"*curr_constraint*"==0.0"

    MOI.set(optimizer, MOI.ConstraintName(), c, curr_constraint)
    append!(constraint_dict["zoneoutflow_c"], [curr_constraint])
end

"""
Add aggregation constraints which aggregate all flows on an arc 
INPUT:
optimizer: optimization model
x : variable
x_agg: aggegation variable
edge: edge being considered
constraint_dict: dictionary of all constraints
g1: graph data strcuture
params: other graph parameters
OUTPUT:
None
"""
function add_aggregation_constraint(optimizer, x, x_agg, edge, constraint_dict, g1, params)
    curr_edge = params.edge_list[edge]

    noz = params.no_zones
    agg_term = [MOI.ScalarAffineTerm(1.0, x_agg[edge])]
    split_terms = [MOI.ScalarAffineTerm(-1.0, x[vl(dest_zone, params.edge_list[edge][1], params.edge_list[edge][2], params)]) for dest_zone in range(1, step=1, length=params.no_zones)]

    agg_term = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(agg_term)
    split_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(split_terms)

    curr_constraint = ""
    
    c = MOI.add_constraint(
            optimizer, 
            MOI.ScalarAffineFunction([agg_term; split_terms], 0.0),
            MOI.EqualTo(0.0)
            )

    for dest_zone in 1:noz
        if curr_constraint == ""
            curr_constraint = curr_constraint*"=x_$(dest_zone),$(curr_edge[1])a$(curr_edge[2])"
        else
            curr_constraint = curr_constraint*"+x_$(dest_zone),$(curr_edge[1])a$(curr_edge[2])"
        end
    end


    curr_constraint = "ax_agg_$(curr_edge[1]) a $(curr_edge[2]) "* curr_constraint

    #MOI.set(optimizer, MOI.ConstraintName(), c, curr_constraint)

    append!(constraint_dict["aggregation"], [curr_constraint])
end


