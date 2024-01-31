struct NetworkLMO <: FrankWolfe.LinearMinimizationOracle
    graph::Graphs.SimpleDiGraph{Int64}
    params::GraphParams
end

function FrankWolfe.compute_extreme_point(
    lmo::NetworkLMO,
    direction
)
    v = all_or_nothing(direction, lmo.params, lmo.graph)
    return v
end