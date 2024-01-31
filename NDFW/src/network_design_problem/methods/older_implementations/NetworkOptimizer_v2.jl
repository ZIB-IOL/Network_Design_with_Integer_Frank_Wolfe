# mutable struct Optimizer <: MOI.AbstractOptimizer
    
#     inner::SCIP.SCIPData
#     reference::SCIP.PtrMap
#     constypes::SCIP.ConsTypeMap
#     binbounds::Dict()
#     start::Dict()

#     function Optimizer(; kwargs...)
#         scip = Ref{Ptr{SCIP.SCIP_}}(C_NULL)
#         SCIP.@SCIP_CALL SCIP.SCIPcreate(scip)
#         @assert scip[] != C_NULL
#         SCIP.@SCIP_CALL SCIP.SCIPincludeDefaultPlugins(scip[])
#         SCIP.@SCIP_CALL SCIP.SCIPcreateProbBasic(scip[], "")

#         return new(
#             SCIP.SCIPData(scip, Dict(), Dict(), 0, 0, Dict(), Dict(), Dict(), Dict(), []),
#             PtrMap(),
#             Graphs.SimpleDiGraph{Int64}(), 
#             GraphParams()
#             )
#     end
# end

mutable struct Optimizer <: MOI.AbstractOptimizer
    true_graph::Graphs.SimpleDiGraph{Int64}
    graph::Graphs.SimpleDiGraph{Int64}
    gparams::GraphParams
    sopt::MOIU.Model
    curr_solution::Vector{Float64}
    curr_binary_soln::Vector{Float64}
    curr_status::String
    integer_variables::Vector{Int64}
    f::Function 
    type_feas_cuts::Vector{String}
    method::String

    function Optimizer(; kwargs...)
        o = new(Graphs.SimpleDiGraph{Int64}(),
                Graphs.SimpleDiGraph{Int64}(), 
                GraphParams(), 
                MOIU.Model{Float64}(),
                [],
                [],
                "Unsolved",
                [],
                x -> x,
                [],
                ""
                )

        return o
    end
end

#FrankWolfe.MathOptLMO(o::Optimizer) = FrankWolfe.MathOptLMO(o.sopt)

free_scip(o::Optimizer) = SCIP.free_scip(o.sopt.inner)

MOI.is_empty(o::Optimizer) = MOI.is_empty(o.sopt)

MOI.empty!(o::Optimizer) = MOI.empty!(o.sopt)



# function Base.show(io::IO, model::Optimizer)
#     return print(io, "NDFW solve with the graph $(model.ptr)")
# end

MOI.get(::Optimizer, ::MOI.SolverName) = "SCIP"

MOI.get(o::Optimizer, param::MOI.RawOptimizerAttribute) = MOI.get(o.sopt, param)

MOI.set(o::Optimizer, param::MOI.RawOptimizerAttribute, value) = MOI.set(o.sopt, param, value)

MOI.supports(o::Optimizer, ::MOI.Silent) = true

MOI.get(o::Optimizer, ::MOI.Silent) = MOI.get(o.sopt, MOI.Silent())

MOI.set(o::Optimizer, ::MOI.Silent, value) = MOI.set(o.sopt, MOI.Silent(), value)

MOI.supports(o::Optimizer, ::MOI.TimeLimitSec) = true

MOI.get(o::Optimizer, ::MOI.TimeLimitSec) = MOI.get(o.sopt, MOI.TimeLimitSec())

MOI.set(o::Optimizer, ::MOI.TimeLimitSec, value) = MOI.set(o.sopt, MOI.TimeLimitSec(), value)

MOI.get(o::Optimizer, ::MOI.NumberOfConstraints{F,S}) where {F,S} = MOI.get(o.sopt, MOI.NumberOfConstraints{F,S}()) 

MOI.get(o::Optimizer, ::MOI.ListOfConstraintTypesPresent) = MOI.get(o.sopt, MOI.ListOfConstraintTypesPresent())

MOI.get(o::Optimizer, ::MOI.ListOfConstraintIndices{F, S}) where {F, S} = MOI.get(o.sopt, MOI.ListOfConstraintIndices{F, S}()) 


function MOI.optimize!(o::Optimizer)
    t1 = MOI.get(o, MOI.ObjectiveFunction{SAF}())
    obj_coeffs = [t1.terms[i].coefficient for i in eachindex(t1.terms)]

    total_scenarios = o.gparams.total_scenarios
    edge_count = Graphs.ne(o.graph)
    dest_var_cnt = edge_count * o.gparams.od_pair_count
    bin_var_cnt = length(o.gparams.removed_edges)

    var_count = total_scenarios * (edge_count + dest_var_cnt) + bin_var_cnt

    if length(obj_coeffs) == 0
        #costs = ones(var_count)
        costs = collect(1.0:var_count)
    else     
        @assert length(obj_coeffs) == var_count   
        costs = obj_coeffs
    end

    if o.method == "NLMO-IFW"
        soln = optimize_over_network(o, costs)
    elseif o.method == "NLMO-P"
        soln = optimize_with_penalty(o, costs)
    else
        throw(ArgumentError("Invalid method"))
    end


    return nothing
end

# "Go back from solved stage to problem modification stage, invalidating results."
allow_modification(o::Optimizer) = SCIP.allow_modification(o.sopt)
