# indices
const VI = MOI.VariableIndex
const CI = MOI.ConstraintIndex
# supported functions
const SAF = MOI.ScalarAffineFunction{Float64}
const SQF = MOI.ScalarQuadraticFunction{Float64}
const VAF = MOI.VectorAffineFunction{Float64}
const VECTOR = MOI.VectorOfVariables
# supported sets
const BOUNDS = Union{MOI.EqualTo{Float64},MOI.GreaterThan{Float64},
    MOI.LessThan{Float64},MOI.Interval{Float64}}
const VAR_TYPES = Union{MOI.ZeroOne,MOI.Integer}
const SOS1 = MOI.SOS1{Float64}
const SOS2 = MOI.SOS2{Float64}
# other MOI types
const AFF_TERM = MOI.ScalarAffineTerm{Float64}
const QUAD_TERM = MOI.ScalarQuadraticTerm{Float64}
const VEC_TERM = MOI.VectorAffineTerm{Float64}