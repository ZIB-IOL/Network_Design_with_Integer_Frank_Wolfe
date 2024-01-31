# Objective function & sense.
#
# SCIP only supports affine objectives. For quadratic or nonlinear objectives,
# the solver will depend on bridges with auxiliary variables. Single variable
# objectives are also accepted, but the type is not correctly remembered.

MOI.supports(::Optimizer, ::MOI.ObjectiveSense) = true
MOI.supports(::Optimizer, ::MOI.ObjectiveFunction{SAF}) = true

MOI.set(o::Optimizer, t::MOI.ObjectiveFunction{SAF}, obj::SAF) = MOI.set(o.sopt, t, obj)

MOI.get(o::Optimizer, t::MOI.ObjectiveFunction{SAF}) = MOI.get(o.sopt, t)

MOI.set(o::Optimizer, t::MOI.ObjectiveSense, sense::MOI.OptimizationSense) = MOI.set(o.sopt, t, sense)

MOI.get(o::Optimizer, s::MOI.ObjectiveSense) = MOI.get(o.sopt, s)

MOI.modify(o::Optimizer, t::MOI.ObjectiveFunction{SAF}, change::MOI.ScalarCoefficientChange{Float64}) = MOI.modify(o.sopt, t, change)

MOI.get(::Optimizer, ::MOI.ObjectiveFunctionType) = SAF
