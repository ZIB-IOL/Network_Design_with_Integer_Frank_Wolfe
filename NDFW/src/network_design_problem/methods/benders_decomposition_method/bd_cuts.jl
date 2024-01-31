"""
This function adds optimality constraints to the bender's master problem given the 
solutions values.
"""
function generate_cut(scen, gparams, y, r, s, t)
    #@info r
    for i in 1:gparams.od_pair_count
        for j in 1:gparams.od_pair_count
            #@info "r[$i][$j] = $(r[i,j])"
            @assert r[i, j] != Inf
        end
    end

    function edge_idx(e)
        edge = gparams.removed_edges[e]
        return gparams.link_dic[edge[1], edge[2]]
    end

    td_2_dest = [maximum([sum(gparams.travel_demand[s][:,z]) for s in 1:gparams.total_scenarios]) for z in 1:gparams.od_pair_count]
    #td_2_dest = [sum(gparams.travel_demand[scen][:, z]) for z in 1:gparams.od_pair_count] #Total travel demand to each destination

    terms = [MOI.ScalarAffineTerm(td_2_dest[z] * s[e, z], y[e]) for e in eachindex(y) for z in 1:gparams.od_pair_count]
    terms_vals = sum([td_2_dest[z] * s[e, z] for e in eachindex(y) for z in 1:gparams.od_pair_count])
    var_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(terms)
    coeff = sum(
        [
        ((r[i, z] - t[z]) * gparams.travel_demand[scen][i, z]) for i in 1:gparams.od_pair_count
        for z in 1:gparams.od_pair_count if gparams.travel_demand[scen][i, z] > 0.0
    ]
    )

    #normalizing_consts = append!(normalizing_consts, coeff - terms_vals)
    #normalizing_consts = append!(normalizing_consts, 0.0)
    #@info "coeff = $coeff"
    #println("r: $(r)")
    #@info "normalizing_consts = $normalizing_consts"

    return var_terms, coeff
end
#TODO allow MAX_FLOW to be different for each destination

"""
This is a simplified version of the feasibility cut. Can be expanded. 
"""
function get_feasibility_cuts(y_vals, y)
    one_sums = 0.0
    terms = []
    for yi in eachindex(y)
        if y_vals[yi] == 0.0
            term = MOI.ScalarAffineTerm(1.0, y[yi])
        elseif y_vals[yi] == 1.0
            term = MOI.ScalarAffineTerm(-1.0, y[yi])
            one_sums += 1.0
        else
            throw(ArgumentError("Incorrect value for y."))
        end
        push!(terms, term)
    end
    #terms = MOI.ScalarAffineTerm.(1, y)
    var_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(terms)
    #coeff = prev_coeff + 1
    return var_terms, 1.0 - one_sums
end

"""
This is a feasibility cut for the subproblem. 
It ensures that at least one edge exists between the origin component 
and the destination component. 
bridge_edges: at least one of these edges must be added to the graph. 
"""
function get_bridge_feasibility_cuts(gparams, y, bridge_edges)
    coeff = 1.0
    terms = []

    #find the index of the bridge edge in the y vector
    for edge in bridge_edges

        yi = findall(x -> x == edge, gparams.removed_edges)[1]

        term = MOI.ScalarAffineTerm(1.0, y[yi])

        push!(terms, term)
    end
    #terms = MOI.ScalarAffineTerm.(1, y)
    var_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(terms)
    #coeff = prev_coeff + 1
    return var_terms, coeff
end

"""
Add optimality cuts to the master problem. 
INPUT
o: Boscia object
bopt: master problem
yl: list of binary variables
η: variable corresponding to the upper bound
r: dual variable corresponding to potential
s: dual variable corresponding to flow
"""
function add_optimality_cuts(scen, gparams, yl, curr_eta, r, s, t)
    opt_cut_var, opt_cut_coeffs = generate_cut(scen, gparams, yl, r, s, t)

    old_opt_coeff = opt_cut_coeffs
    eta_term = MOI.ScalarAffineTerm(1.0, curr_eta)
    push!(opt_cut_var, eta_term)

    #@info "Adding cut: $(opt_cut_coeffs) <= sum of y_i)"
    # MOI.add_constraint(
    #     bopt,
    #     MOI.ScalarAffineFunction(opt_cut_var, 0.0),
    #     MOI.GreaterThan(opt_cut_coeffs)
    # )

    c = (MOI.ScalarAffineFunction(opt_cut_var, 0.0), MOI.GreaterThan(opt_cut_coeffs))

    return [c]
end




"""
This function adds the feasibility cuts to the model. 
INPUT 
o: Boscia object
yl: list of binary variables corresponding to removed edges
state: state of the graph
costs: costs of the edges
"""
function add_feasibility_cuts(yl, yvals, state, src, costs, orig_edge_list, gparams, graph, type_feas_cuts)
    feasibility_cuts = []
    if "bridge" in type_feas_cuts
        bridge_edges = get_bridge_edges(gparams, state, src, costs, orig_edge_list)

        feas_terms, feas_coeff = get_bridge_feasibility_cuts(gparams, yl, bridge_edges)

        # c = MOI.add_constraint(
        #     bopt,
        #     MOI.ScalarAffineFunction(feas_terms, 0.0),
        #     MOI.GreaterThan(feas_coeff)
        # )
        c = (MOI.ScalarAffineFunction(feas_terms, 0.0), MOI.GreaterThan(feas_coeff))

        push!(feasibility_cuts, c)
    end

    if "multi" in type_feas_cuts
        bridge_edges = get_multi_step_bridge_edges(gparams, graph, state, src, costs, orig_edge_list)

        for step in keys(bridge_edges)
            be = bridge_edges[step]
            feas_terms, feas_coeff = get_bridge_feasibility_cuts(gparams, yl, be)

            # c = MOI.add_constraint(
            #     bopt,
            #     MOI.ScalarAffineFunction(feas_terms, 0.0),
            #     MOI.GreaterThan(feas_coeff)
            # )

            c = (MOI.ScalarAffineFunction(feas_terms, 0.0), MOI.GreaterThan(feas_coeff))

            push!(feasibility_cuts, c)
        end
    end

    if "point" in type_feas_cuts
        feas_terms, feas_coeff = get_feasibility_cuts(yvals, yl)

        # c = MOI.add_constraint(
        #     bopt,
        #     MOI.ScalarAffineFunction(feas_terms, 0.0),
        #     MOI.GreaterThan(feas_coeff)
        # )

        c = (MOI.ScalarAffineFunction(feas_terms, 0.0), MOI.GreaterThan(feas_coeff))

        push!(feasibility_cuts, c)
    end

    return feasibility_cuts
end


function add_feasibility_cuts(gparams, yl, r, t, s, scen)
    feasibility_cuts = []

    td_2_dest = [maximum([sum(gparams.travel_demand[s][:, z]) for s in 1:gparams.total_scenarios]) for z in 1:gparams.od_pair_count]
    #td_2_dest = [sum(gparams.travel_demand[scen][:, z]) for z in 1:gparams.od_pair_count] #Total travel demand to each destination

    terms = [MOI.ScalarAffineTerm(td_2_dest[z] * s[e, z], yl[e]) for e in eachindex(yl) for z in 1:gparams.od_pair_count]
    feas_terms = Vector{MathOptInterface.ScalarAffineTerm{Float64}}(terms)

    feas_coeff = sum(
        [
        ((r[i, z] - t[z]) * gparams.travel_demand[scen][i, z]) for i in 1:gparams.od_pair_count
        for z in 1:gparams.od_pair_count if gparams.travel_demand[scen][i, z] > 0.0
    ]
    )

    # c = MOI.add_constraint(
    #     bopt,
    #     MOI.ScalarAffineFunction(feas_terms, 0.0),
    #     MOI.GreaterThan(feas_coeff)
    # )

    c = (MOI.ScalarAffineFunction(feas_terms, 0.0), MOI.GreaterThan(feas_coeff))

    push!(feasibility_cuts, c)

    return feasibility_cuts
end
