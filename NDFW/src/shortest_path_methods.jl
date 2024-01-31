
function TA_dijkstra_shortest_paths_default(graph, travel_time, origin, init_node, term_node, first_thru_node, distmx)
    no_arc = ne(graph)

    distmx .= Inf
    for i in 1:no_arc
      if term_node[i] >= first_thru_node
          distmx[init_node[i], term_node[i]] = travel_time[i]
      end
    end

    state = dijkstra_shortest_paths(graph, origin, distmx)
    return state
end

function TA_reverse_shortest_paths(graph, travel_time, origin, link_dic, first_thru_node, distmx)
    #no_node = Graphs.nv(graph)
    #no_arc = Graphs.ne(graph)

    max_factor = 100 * sum(travel_time)
    distmx .= max_factor
    curr_edge_list = collect(Graphs.edges(graph))

    for edge in curr_edge_list
        src = edge.src
        dst = edge.dst

        #TODO get better distance matrix for shortest path
        if src >= first_thru_node 
            distmx[src, dst] = travel_time[link_dic[src, dst]]
        end
    end

    # for i in 1:no_arc
    #     if term_node[i] >= first_thru_node
    #         distmx[init_node[i], term_node[i]] = travel_time[i]
    #     end
    # end

    state = Graphs.dijkstra_shortest_paths(graph, origin, distmx)

    #Update the distances and parents to reflect the actual travel times
    #This means travelling through at least 2 arcs with time exceed the max_factor

    state.parents[state.dists.>2*max_factor] .= 0
    state.dists[state.dists.>2*max_factor] .= Inf

    return state
end


function TA_dijkstra_shortest_paths(graph, travel_time, origin, link_dic, first_thru_node, distmx)
    #no_node = Graphs.nv(graph)
    #no_arc = Graphs.ne(graph)

    max_factor = 100 * sum(travel_time)
    distmx .= max_factor

    curr_edge_list = collect(Graphs.edges(graph))
    
    for edge in curr_edge_list
        src = edge.src
        dst = edge.dst

        #TODO get better distance matrix for shortest path
        if dst >= first_thru_node
            distmx[src, dst] = travel_time[link_dic[src, dst]]
        end

        if dst < first_thru_node
            distmx[src, dst] = max_factor + travel_time[link_dic[src, dst]]
        end
    end
    
    # for i in 1:no_arc
    #     if term_node[i] >= first_thru_node
    #         distmx[init_node[i], term_node[i]] = travel_time[i]
    #     end
    # end

    state = Graphs.dijkstra_shortest_paths(graph, origin, distmx)

    state.parents[state.dists.>2*max_factor] .= 0
    state.dists[state.dists.>2*max_factor] .= Inf

    return state
end

function TA_floyd_warshall_shortest_paths(graph, travel_time, link_dic, first_thru_node)
    no_node = Graphs.nv(graph)
    no_arc = Graphs.ne(graph)

    distmx = 1e6*ones(no_node, no_node)
    curr_edge_list = collect(Graphs.edges(graph))
    
    for edge in curr_edge_list
        src = edge.src
        dst = edge.dst

        #TODO get better distance matrix for shortest path
        if dst >= first_thru_node
            distmx[src, dst] = travel_time[link_dic[src, dst]]
        elseif dst < first_thru_node
            distmx[src, dst] = 1e6
        end
    end
    
    # for i in 1:no_arc
    #     if term_node[i] >= first_thru_node
    #         distmx[init_node[i], term_node[i]] = travel_time[i]
    #     end
    # end

    state = Graphs.floyd_warshall_shortest_paths(graph, distmx)
    return state
end

function TA_dijkstra_shortest_paths(graph, travel_time, origin, link_dic)
    no_node = nv(graph)
    no_arc = ne(graph)

    distmx = Inf*ones(no_node, no_node)

    curr_edge_list = collect(Graphs.edges(graph))
    
    for edge in curr_edge_list
        src = edge.src
        dst = edge.dst
        distmx[src, dst] = travel_time[link_dic[src, dst]]
    end

    # for i in 1:no_arc
    #   distmx[init_node[i], term_node[i]] = travel_time[i]
    # end

    state = Graphs.dijkstra_shortest_paths(graph, origin, distmx)
    return state
end

function create_graph(init_node, term_node)
    @assert Base.length(init_node)==Base.length(term_node)

    no_node = max(maximum(init_node), maximum(term_node))
    no_arc = Base.length(init_node)

    graph = Graphs.SimpleDiGraph(no_node)
    for i=1:no_arc
        Graphs.add_edge!(graph, init_node[i], term_node[i])
    end
    return graph
end



function get_vector(state, origin, destination, link_dic)
    current = destination
    parent = -1
    x = zeros(Int, maximum(link_dic))

    while parent != origin && origin != destination && current != 0
        parent = state.parents[current]

        # println("origin=$origin, destination=$destination, parent=$parent, current=$current")

        if parent != 0
            link_idx = link_dic[parent,current]
            if link_idx != 0
                x[link_idx] = 1
            end
        end

        current = parent
    end

    return x
end

function add_demand_vector!(x, demand, state, origin, destination, link_dic)
  current = destination
  parent = -1

  while parent != origin && origin != destination && current != 0
      parent = state.parents[current]

      if parent != 0
          link_idx = link_dic[parent,current]
          if link_idx > length(x)
            println("parent: $parent, current: $current, link_idx: $link_idx")
          end
          if link_idx != 0
              x[link_idx] += demand
          end
      end

      current = parent
  end
end

function add_demand_vector_detailed!(
    x::Vector{Float64}, 
    demand::Float64, 
    state::Graphs.DijkstraState{Float64, Int64}, 
    origin::Int64, 
    destination::Int64, 
    link_dic::SparseMatrixCSC{Int64, Int64}, 
    edge_list::Vector{Tuple{Int64, Int64}}, 
    od_pair_count::Int64
    )

    current = destination
    parent = -1
  
    agg_arc_start = length(edge_list) * od_pair_count

    while parent != origin && origin != destination && current != 0
        parent = state.parents[current]
  
        if parent != 0
            link_idx = link_dic[parent,current]
            if link_idx != 0
                x[(destination - 1) * length(edge_list) + link_idx] += demand
                x[agg_arc_start + link_idx] += demand
            end
        end
  
        current = parent
    end
end

function add_demand_vector_detailed_reverse!(
    x::Vector{Float64},
    demand::Float64,
    state::Graphs.DijkstraState{Float64,Int64},
    origin::Int64,
    destination::Int64,
    link_dic::SparseMatrixCSC{Int64,Int64},
    edge_list::Vector{Tuple{Int64,Int64}},
    od_pair_count::Int64
)

    current = destination
    parent = -1

    agg_arc_start = length(edge_list) * od_pair_count

    while parent != origin && origin != destination && current != 0
        parent = state.parents[current]

        if parent != 0
            link_idx = link_dic[parent, current]
            if link_idx != 0
                x[(origin-1)*length(edge_list)+link_idx] += demand
                x[agg_arc_start+link_idx] += demand
            end
        end

        current = parent
    end
end


function get_path(state, origin, destination)
    current = destination
    parent = -1
    path = []

    while parent != origin && origin != destination && current != 0
        parent = state.parents[current]
        push!(path, (parent, current))
        current = parent
    end

    return [reverse(path)]
end

function get_all_shortest_paths(state, origin, destination)
    
    all_paths = []
    invalid_paths = []
    valid_paths = state.pathcounts[destination]
    while length(all_paths) + length(invalid_paths) < state.pathcounts[destination]
        @info "current path count: $(length(all_paths))"
        current = destination
        pred = -1
        path = []

        temp_path = [] #to check for loops
        
        current_path_valid = true
        while pred != origin && origin != destination && current != 0 && pred != destination
            println("current: $current, pred: $pred, origin: $origin, destination: $destination")
            println("temp path $(temp_path)")
            preds = state.predecessors[current]
            if length(preds) == 1
                pred = preds[1]
                push!(path, (pred, current))
                push!(temp_path, (pred, current))
                current = pred

                if pred == destination
                    @info "looping back to destination"
                    #valid_paths -= 1
                    current_path_valid = false
                    break
                elseif (pred, current) in temp_path
                    break
                end

            elseif length(preds) > 1
                for j in eachindex(preds)
                    pred = preds[j]
                    if pred == destination
                        @info "looping back to destination"
                        #valid_paths -= 1
                        current = 0
                        current_path_valid = false
                        break
                    elseif (pred, current) in temp_path
                        current = 0
                        break
                    end

                    if edge_in_path((pred, current), all_paths) || edge_in_path((pred, current), invalid_paths)
                        continue
                    else
                        push!(path, (pred, current))
                        push!(temp_path, (pred, current))
                        current = pred
                        break
                    end
                end
            else
                throw(
                    ArgumentError(
                        "No predecessor found for node $current"
                    )
                )
            end
        end

        if current_path_valid
            push!(all_paths, reverse(path))
        else
            push!(invalid_paths, reverse(path))
        end

        if length(all_paths) >= valid_paths
            break
        end
    end

    return all_paths, invalid_paths
end


function get_all_shortest_paths2(state, origin, destination)
    
    all_paths = []
    invalid_paths = []

    path_seeds = Vector{Any}([[destination]])

    while length(path_seeds) > 0
        seed = pop!(path_seeds)

    
        current = seed[end]
        seed_path = seed

        pred = -1
        path = []



        while pred != origin && origin != destination && current != 0 && pred != destination && !(pred in seed_path[1:(end - 1)])
            preds = state.predecessors[current]

            if length(preds) == 1
                pred = preds[1]
                #push!(path, (pred, current))
                push!(seed_path, pred)
                path = [(seed_path[i],seed_path[i+1]) for i in 1:length(seed_path)-1]
                current = pred
            elseif length(preds) > 1
                for pred in preds
                    #push!(seed_path, pred)
                    push!(path_seeds, [seed_path; pred])
                end
                break
            end
        end

        if seed_path[end] == origin
            push!(all_paths, reverse(path))
        end

        paths_to_delete = []
        for pi in eachindex(path_seeds)
            p = path_seeds[pi]
            if p[end] == destination
                push!(invalid_paths, reverse(p))
                append!(paths_to_delete, pi)
            elseif p[end] in p[1:(end - 1)]
                push!(invalid_paths, reverse(p))
                append!(paths_to_delete, pi)
            end
        end

        deleteat!(path_seeds, paths_to_delete)
        

    end

    return all_paths, invalid_paths
end

function edge_in_path(edge, paths)
    for path in paths
        if edge in path
            return true
        end
    end
    return false
end






# function get_shortest_path(init_node, term_node, link_length, origin, destination)
#     @assert Base.length(init_node)==Base.length(term_node)
#     @assert Base.length(init_node)==Base.length(link_length)
#
#     graph = create_graph(init_node, term_node)
#
#     state = dijkstra_shortest_paths(graph, link_length, origin)
#
#     path = get_path(state, origin, destination)
#     x = get_vector(path, init_node, term_node)
#
#     return path, x
# end
#
# function get_vector(state, origin, destination, init_node, term_node)
#     current = destination
#     parent = -1
#     x = zeros(Int, Base.length(init_node))
#
#     while parent != origin
#         parent = state.parents[current]
#
#         for j=1:Base.length(init_node)
#             if init_node[j]==parent && term_node[j]==current
#                 x[j] = 1
#                 break
#             end
#         end
#
#         current = parent
#     end
#
#     return x
# end

