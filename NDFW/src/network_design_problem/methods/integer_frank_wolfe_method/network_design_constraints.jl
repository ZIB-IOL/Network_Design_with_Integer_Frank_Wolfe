"""
Add network design constraints dependigng upon specified arcs. 
INPUT:
optimizer: optimization model
new_edge_list: edges on which to add network design constraints
constraint_dict: dictionary of all constraints
g1: graph data strcuture
params: other graph parameters
OUTPUT:
None
"""
function add_network_design_constraints!(
    optimizer,
    new_edge_list,
    constraint_dict,
    g1::MC_graph_with_weights, 
    params::Params,
    total_scenarios::Int64,
    variables_initialized = false,
    y = nothing
)
    new_edge_count = length(new_edge_list)

    constraint_dict["network_design"] = []

    #check if binary variables are initialized or provided
    if variables_initialized == false
        y = initialize_nd_variables(optimizer, new_edge_count, new_edge_list)
    else
        @assert isnothing(y) == false "y is not initialized"
    end


    td_2_dest = [maximum([sum(params.node_demands_dict[s][z]) for s in 1:total_scenarios]) for z in 1:params.no_zones]
    #add network design constraints for every scenario
    for scen in 1:total_scenarios

        var_cnt = params.no_edges * params.no_zones
        total_scen_vars = params.no_edges * params.no_zones + params.no_edges

        first_var = (scen-1) * total_scen_vars

        x = [MOI.VariableIndex(i) for i=(first_var + 1) : (first_var + var_cnt)]
        #x_agg = [MOI.VariableIndex(i) for i=(first_var + var_cnt+1):(first_var + var_cnt + params.no_edges)]
        
        

        for (ec, curr_edge) in enumerate(new_edge_list)
            for dest_zone in 1:params.no_zones
                src, dst = curr_edge
                edge = vl(1, src, dst, params)
                split_term = MOI.ScalarAffineTerm(1.0, x[vl(dest_zone, params.edge_list[edge][1], params.edge_list[edge][2], params)]) 

                flow_term = [split_term]
                bound_term = [MOI.ScalarAffineTerm(-td_2_dest[dest_zone], y[ec])]

                MOI.add_constraint(
                    optimizer, 
                    MOI.ScalarAffineFunction([flow_term; bound_term], 0.0),
                    MOI.LessThan(0.0)
                    )

                curr_constraint = ""

                curr_constraint = curr_constraint*"x_$(dest_zone),$(curr_edge[1]) a $(curr_edge[2])"

                curr_constraint = curr_constraint*" - $(params.MAX_FLOW).y_$(curr_edge[1]) a $(curr_edge[2]) <= 0" 

                append!(constraint_dict["network_design"], [curr_constraint])
            end
        end

        #@info "Added $(length(constraint_dict["network_design"])) network design constraints. "
    end

    if variables_initialized == false
        return y
    end
end



"""
    add_network_design_constraints_inverted!(
        optimizer,
        new_edge_list,
        constraint_dict,
        g1::MC_graph_with_weights,
        params::Params,
        total_scenarios::Int64,
        variables_initialized=false,
        y=nothing
    )

Add network design constraints for the inverted integer Frank-Wolfe method.
This function uses bigM constraints to impose the network design constraints.
However it imposes them in the form x <= M (1 - y) where y is the binary variable. 
This is necessary for comparison against the indicator constraints.

# Arguments
- `optimizer`: The optimizer object.
- `new_edge_list`: A list of new edges.
- `constraint_dict`: A dictionary to store the constraints.
- `g1::MC_graph_with_weights`: The graph object.
- `params::Params`: The parameters object.
- `total_scenarios::Int64`: The total number of scenarios.
- `variables_initialized`: A boolean indicating whether the binary variables are initialized or provided. Default is `false`.
- `y`: The binary variables. Default is `nothing`.

# Returns
- If `variables_initialized` is `false`, returns the binary variables `y`.

"""
function add_network_design_constraints_inverted!(
    optimizer,
    new_edge_list,
    constraint_dict,
    g1::MC_graph_with_weights,
    params::Params,
    total_scenarios::Int64,
    variables_initialized=false,
    y=nothing
)
    new_edge_count = length(new_edge_list)

    constraint_dict["network_design"] = []

    #check if binary variables are initialized or provided
    if variables_initialized == false
        y = initialize_nd_variables(optimizer, new_edge_count, new_edge_list)
    else
        @assert isnothing(y) == false "y is not initialized"
    end


    td_2_dest = [maximum([sum(params.node_demands_dict[s][z]) for s in 1:total_scenarios]) for z in 1:params.no_zones]
    #add network design constraints for every scenario
    for scen in 1:total_scenarios

        var_cnt = params.no_edges * params.no_zones
        total_scen_vars = params.no_edges * params.no_zones + params.no_edges

        first_var = (scen - 1) * total_scen_vars

        x = [MOI.VariableIndex(i) for i = (first_var+1):(first_var+var_cnt)]
        #x_agg = [MOI.VariableIndex(i) for i=(first_var + var_cnt+1):(first_var + var_cnt + params.no_edges)]



        for (ec, curr_edge) in enumerate(new_edge_list)
            for dest_zone in 1:params.no_zones
                src, dst = curr_edge
                edge = vl(1, src, dst, params)
                split_term = MOI.ScalarAffineTerm(1.0, x[vl(dest_zone, params.edge_list[edge][1], params.edge_list[edge][2], params)])

                flow_term = [split_term]
                bound_term = [MOI.ScalarAffineTerm(-td_2_dest[dest_zone], y[ec])]

                MOI.add_constraint(
                    optimizer,
                    MOI.ScalarAffineFunction([flow_term; bound_term], 0.0),
                    MOI.LessThan(td_2_dest[dest_zone])
                )

                curr_constraint = ""

                curr_constraint = curr_constraint * "x_$(dest_zone),$(curr_edge[1]) a $(curr_edge[2])"

                curr_constraint = curr_constraint * " - $(params.MAX_FLOW).y_$(curr_edge[1]) a $(curr_edge[2]) <= 0"

                append!(constraint_dict["network_design"], [curr_constraint])
            end
        end

        #@info "Added $(length(constraint_dict["network_design"])) network design constraints. "
    end

    if variables_initialized == false
        return y
    end
end



"""
    add_network_design_constraints_indicator!(
        optimizer,
        new_edge_list,
        constraint_dict,
        g1::MC_graph_with_weights,
        params::Params,
        total_scenarios::Int64,
        variables_initialized = false,
        y = nothing
    )

Add network design constraints for the integer Frank-Wolfe method.
This particular function uses indicator constraints to impose the network design constraints. 
As of now this only works with SCIP. Gurobi throws bugs


# Arguments
- `optimizer`: The optimizer object.
- `new_edge_list`: A list of new edges.
- `constraint_dict`: A dictionary to store the constraints.
- `g1::MC_graph_with_weights`: The graph object.
- `params::Params`: The parameters object.
- `total_scenarios::Int64`: The total number of scenarios.
- `variables_initialized`: A boolean indicating whether the binary variables are initialized or provided. Default is `false`.
- `y`: The binary variables. Default is `nothing`.

# Returns
- If `variables_initialized` is `false`, returns the initialized binary variables `y`.

# Example
"""
function add_network_design_constraints_indicator!(
    optimizer,
    new_edge_list,
    constraint_dict,
    g1::MC_graph_with_weights,
    params::Params,
    total_scenarios::Int64,
    variables_initialized = false,
    y = nothing
)


    new_edge_count = length(new_edge_list)

    constraint_dict["network_design"] = []

    #check if binary variables are initialized or provided
    if variables_initialized == false
        y = initialize_nd_variables(optimizer, new_edge_count, new_edge_list)
    else
        @assert isnothing(y) == false "y is not initialized"
    end


    #add network design constraints for every scenario
    for scen in 1:total_scenarios

        var_cnt = params.no_edges * params.no_zones
        total_scen_vars = params.no_edges * params.no_zones + params.no_edges
        #x_agg = [MOI.VariableIndex(i) for i=(var_cnt+1):(var_cnt + params.no_edges)]
        
        first_var = (scen - 1) * total_scen_vars

        x = [MOI.VariableIndex(i) for i = (first_var+1):(first_var+var_cnt)]
        #x_agg = [MOI.VariableIndex(i) for i = (first_var+var_cnt+1):(first_var+var_cnt+params.no_edges)]

        #@info "Loading aggrgate edge variables $(x_agg[1]) to $(x_agg[end])"

        for (ec, curr_edge) in enumerate(new_edge_list)
            for dest_zone in 1:params.no_zones
                src, dst = curr_edge
                edge = vl(1, src, dst, params)
                split_term = MOI.ScalarAffineTerm(1.0, x[vl(dest_zone, params.edge_list[edge][1], params.edge_list[edge][2], params)])

                #flow_term = [MOI.ScalarAffineTerm(1.0, x_agg[edge])]
                #bound_term = [MOI.ScalarAffineTerm(-params.MAX_FLOW, y[ec])]

                gl = MOI.VectorAffineFunction(
                    [
                        MOI.VectorAffineTerm(1, MOI.ScalarAffineTerm(1.0, y[ec])),
                        MOI.VectorAffineTerm(2, split_term),
                    ],
                    [0.0, 0.0]
                )

                c = MOI.add_constraint(
                    optimizer, 
                    gl,
                    MOI.Indicator{MOI.ACTIVATE_ON_ONE}(MOI.LessThan(0.0))
                    )

                curr_constraint = ""

                curr_constraint = curr_constraint*"x_agg_$(curr_edge[1]) a $(curr_edge[2])"

                curr_constraint = curr_constraint*" - $(params.MAX_FLOW).y_$(curr_edge[1]) a $(curr_edge[2]) <= 0" 

                append!(constraint_dict["network_design"], [curr_constraint])

                #println("added constraint $(c)")
            end
        end
    end

    @info "Added $(length(constraint_dict["network_design"])) network design constraints. "

    if variables_initialized == false
        return y
    end
end


"""
    initialize_nd_variables(optimizer, new_edge_count, new_edge_list)

Initialize network design variables for the integer Frank-Wolfe method.

# Arguments
- `optimizer`: The optimizer object.
- `new_edge_count`: The number of new edges.
- `new_edge_list`: A list of tuples representing the source and destination of each new edge.

# Returns
- `y`: An array of network design variables.

"""
function initialize_nd_variables(optimizer, new_edge_count, new_edge_list)
    y = MOI.add_variables(optimizer, new_edge_count)

    #@info "Initialized $(new_edge_count) network design variables from $(y[1]) to $(y[end])"

    for i in 1:new_edge_count
        MOI.add_constraint(optimizer, y[i], MOI.LessThan(1.0))
        MOI.add_constraint(optimizer, y[i], MOI.GreaterThan(0.0))
        MOI.add_constraint(optimizer, y[i], MOI.ZeroOne())
        zone = 1
        (src, dst) = new_edge_list[i]
        MOI.set(optimizer, MOI.VariableName(), y[i], "y_$(zone),$(src) a $(dst)")
    end
   
    return y
end
