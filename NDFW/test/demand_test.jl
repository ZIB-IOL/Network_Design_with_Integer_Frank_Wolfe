function check_demand_for_node(x, gparams, g)
    incoming, outgoing = NDFW.neighbouring_edges(g)
    edge_count = length(gparams.init_nodes)
    for scen in 1:gparams.total_scenarios
        for dest in 1:gparams.od_pair_count
            soln = sum(x[(dest-1)*edge_count.+incoming[dest]])
            demand = sum(gparams.travel_demand[scen][:, dest])
            @assert soln ≈ demand "$(dest) for $(scen), $(soln) vs $(demand)"
        end
    end

    for scen in 1:gparams.total_scenarios
        for src in 1:gparams.od_pair_count
            for dest in 1:gparams.od_pair_count
                soln = sum(x[(dest-1)*edge_count.+outgoing[src]])
                demand = sum(gparams.travel_demand[scen][src, dest])
                @assert soln ≈ demand "$(src) -> $(dest) for $(scen), $(soln) vs $(demand)"
            end
        end
    end
end

function balance_test(x, gparams, g)
    incoming, outgoing = NDFW.neighbouring_edges(g)
    edge_count = length(gparams.init_nodes)
    nodes = Graphs.nv(g)

    for scen in 1:gparams.total_scenarios
        for dest in 1:gparams.od_pair_count
            for node in (gparams.od_pair_count+1):nodes
                inc = 0.0
                out = 0.0
                if node in keys(incoming)
                    inc += sum(x[(dest-1)*edge_count.+incoming[node]])
                end

                if node in keys(outgoing)
                    out += sum(x[(dest-1)*edge_count.+outgoing[node]])
                end

                @assert inc ≈ out "$(dest) for node $(node) for $(scen), $(inc) vs $(out)"
            end
        end
    end
end

        
    