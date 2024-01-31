## variables (general)


MOI.add_variable(o::Optimizer) = MOI.add_variable(o.sopt)

MOI.add_variables(o::Optimizer, n) = MOI.add_variables(o.sopt, n)
MOI.get(o::Optimizer, nov::MOI.NumberOfVariables) = MOI.get(o.sopt, nov)
MOI.get(o::Optimizer, lvi::MOI.ListOfVariableIndices) = MOI.get(o.sopt, lvi)
MOI.is_valid(o::Optimizer, vi::MOI.VariableIndex) = MOI.is_valid(o.sopt, vi)

MOI.get(o::Optimizer, vn::MOI.VariableName, vi::MOI.VariableIndex)::String = MOI.get(o.sopt, vn, vi)::String

MOI.set(o::Optimizer, vn::MOI.VariableName, vi::MOI.VariableIndex, name::String) = MOI.set(o.sopt, vn, vi, name)

MOI.supports(::Optimizer, vn::MOI.VariableName, ::Type{MOI.VariableIndex}) = true

MOI.get(o::Optimizer, t::Type{MOI.VariableIndex}, name::String) = MOI.get(o.sopt, t, name)

MOI.delete(o::Optimizer, vi::MOI.VariableIndex) = MOI.delete(o.sopt, vi)

## variable types (binary, integer)

MOI.supports_constraint(o::Optimizer, ::Type{MOI.VariableIndex}, ::Type{<:VAR_TYPES}) = true

scip_vartype(::Type{MOI.ZeroOne}) = SCIP.SCIP_VARTYPE_BINARY
scip_vartype(::Type{MOI.Integer}) = SCIP.SCIP_VARTYPE_INTEGER

MOI.add_constraint(o::Optimizer, vi::MOI.VariableIndex, set::S) where {S <: VAR_TYPES} = MOI.add_constraint(o.sopt, vi, set) 

MOI.delete(o::Optimizer, ci::MOI.ConstraintIndex{MOI.VariableIndex,S}) where {S <: VAR_TYPES} = MOI.delete(o.sopt, ci) 

## variable bounds

MOI.supports_constraint(o::Optimizer, ::Type{MOI.VariableIndex}, ::Type{<:BOUNDS}) = true

MOI.add_constraint(o::Optimizer, vi::MOI.VariableIndex, set::S) where {S <: BOUNDS} = MOI.add_constraint(o.sopt, vi, set) 

reset_bounds(o::Optimizer, v, lb, ub, ci::MOI.ConstraintIndex{MOI.VariableIndex, S}) where S <: Union{MOI.Interval{Float64}, MOI.EqualTo{Float64}} = reset_bounds(o.sopt, v, lb, ub, ci)

reset_bounds(o::Optimizer, v, lb, ub, ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.GreaterThan{Float64}}) = reset_bounds(o.sopt, v, lb, ub, ci)

reset_bounds(o::Optimizer, v, lb, ub, ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.LessThan{Float64}}) = reset_bounds(o.sopt, v, lb, ub, ci)

MOI.delete(o::Optimizer, ci::MOI.ConstraintIndex{MOI.VariableIndex,S}) where {S <: BOUNDS} = MOI.delete(o.sopt, ci) 

#MOI.set(o::Optimizer, cs::MOI.ConstraintSet, ci::MOI.ConstraintIndex{MOI.VariableIndex,S}, set::S) where {S <: BOUNDS} = MOI.set(o.sopt, cs, ci, set) 

function MOI.set(o::Optimizer, cs::MOI.ConstraintSet, ci::MOI.ConstraintIndex{MOI.VariableIndex,S}, set::S) where {S <: BOUNDS}
    MOI.set(o.sopt, cs, ci, set) 

    #allow_modification(o)
    lb, ub = SCIP.bounds(set)

    edge_count = length(o.gparams.init_nodes)
    od_pair_count = o.gparams.od_pair_count
    total_scenarios = o.gparams.total_scenarios

    total_scenario_vars = od_pair_count * edge_count + edge_count
    
    if lb == 1.0 
        #println("fixing value of variable $ci.value to upper bound $lb")
        o.curr_binary_soln[ci.value - total_scenario_vars * total_scenarios] = lb
        o.curr_solution[ci.value] = lb
    elseif ub == 0.0
        #println("fixing value of variable $ci.value to lower bound $ub")
        o.curr_binary_soln[ci.value - total_scenario_vars * total_scenarios] = ub
        o.curr_solution[ci.value] = ub
    end
    
    #println("curr binary soln: $(o.curr_binary_soln)")
    #update_graph(o, ci.value, lb, ub)

    return nothing
end

# TODO: is actually wrong for unbounded variables?
MOI.is_valid(o::Optimizer, ci::MOI.ConstraintIndex{MOI.VariableIndex,S}) where {S <: BOUNDS} = MOI.is_valid(o.sopt, ci)

MOI.get(o::Optimizer, cf::MOI.ConstraintFunction, ci::MOI.ConstraintIndex{MOI.VariableIndex, S}) where {S <: BOUNDS} = MOI.get(o.sopt, cf, ci)

MOI.get(o::Optimizer, cs::MOI.ConstraintSet, ci::MOI.ConstraintIndex{MOI.VariableIndex, S}) where {S <: BOUNDS} = MOI.get(o.sopt, cs, ci)

MOI.is_valid(o::Optimizer, ci::MOI.ConstraintIndex{MOI.VariableIndex,S}) where {S <: Union{MOI.ZeroOne, MOI.Integer}} = MOI.is_valid(o.sopt, ci)

MOI.get(o::Optimizer, cs::MOI.ConstraintSet, ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.ZeroOne}) = MOI.get(o.sopt, cs, ci)

MOI.get(o::Optimizer, cs::MOI.ConstraintSet, ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.Integer}) = MOI.get(o.sopt, cs, ci)


MOI.get(o::Optimizer, cf::MOI.ConstraintFunction, ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.ZeroOne}) = MOI.get(o.sopt, cf, ci)


MOI.get(o::Optimizer, cf::MOI.ConstraintFunction, ci::MOI.ConstraintIndex{MOI.VariableIndex, MOI.Integer}) = MOI.get(o.sopt, cf, ci)


# # (partial) warm starts

MOI.supports(::Optimizer, vps::MOI.VariablePrimalStart, ::Type{MOI.VariableIndex}) = true

MOI.get(o::Optimizer, vps::MOI.VariablePrimalStart, vi::MOI.VariableIndex) = MOI.get(o.sopt, vps, vi)

MOI.set(o::Optimizer, vps::MOI.VariablePrimalStart, vi::MOI.VariableIndex, value::Float64) = MOI.set(o.sopt, vps, vi, value)


MOI.set(o::Optimizer, vps::MOI.VariablePrimalStart, vi::MOI.VariableIndex, value::Nothing) = MOI.set(o.sopt, vps, vi, value)


MOI.set(::Optimizer, cf::MOI.ConstraintFunction, ci::MOI.ConstraintIndex{MOI.VariableIndex}, vi::MOI.VariableIndex) = MOI.set(o.sopt, cf, ci, vi)

