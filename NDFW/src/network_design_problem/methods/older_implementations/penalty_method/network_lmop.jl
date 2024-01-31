#LMO for the penalty approach
struct NetworkLMO_P <: FrankWolfe.LinearMinimizationOracle
    o::Optimizer
    graph::Graphs.SimpleDiGraph{Int64}
    params::GraphParams
end

NetworkLMO_P(o::Optimizer) = NetworkLMO_P(o, o.graph, o.gparams)

function FrankWolfe.compute_extreme_point(
    lmo::NetworkLMO_P,
    direction;
    kwargs...
)

    soln = optimize_with_penalty(lmo.o, direction)
    
    return soln
end