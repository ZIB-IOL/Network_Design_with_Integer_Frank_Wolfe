struct SPLMO <: FrankWolfe.LinearMinimizationOracle
    graph
    gparams::GraphParams
    distmx::Matrix{Float64}
end

SPLMO(graph, gparams) = SPLMO(graph, gparams, zeros(Float64, Graphs.nv(graph), Graphs.nv(graph)))


function FrankWolfe.compute_extreme_point(splmo::SPLMO, direction; kwargs...)
    graph = splmo.graph
    gparams = splmo.gparams

    ie = Graphs.ne(graph)
    distmx = splmo.distmx

    costs = direction


    curr_solution = get_flow_extreme_point(costs, graph, ie, gparams, distmx)
    
    return curr_solution
end

function get_flow_extreme_point(costs, graph, ie, gparams, distmx)
    link_dic = gparams.link_dic

    scen = 1
    
    v = all_or_nothing_sp(costs, gparams, graph, scen, link_dic, distmx)
    

    return v
end