struct State
    iter::Int64
    x::Vector{Float64}
    obj::Float64
    dual_gap::Float64
    time::Float64
end



function mod_BPR(x::Vector{Float64}, td::NDFW.TA_Data, norm_const=1)
    bpr = similar(x)
    for i in 1:length(bpr)
        bpr[i] = td.free_flow_time[i] * (1.0 + td.b[i] * (x[i] / td.capacity[i])^td.power[i])
        bpr[i] += td.toll_factor * td.toll[i] + td.distance_factor * td.link_length[i]
    end
    return bpr / norm_const
end

mod_gradient = mod_BPR

function mod_objective(x::Vector{Float64}, td::NDFW.TA_Data, norm_const=1)
    # value = free_flow_time .* ( x + b.* ( x.^(power+1)) ./ (capacity.^power) ./ (power+1))
    # return sum(value)

    sum = 0.0
    for i in 1:length(x)
        sum += td.free_flow_time[i] * (x[i] + td.b[i] * (x[i]^(td.power[i] + 1)) / (td.capacity[i]^td.power[i]) / (td.power[i] + 1))
        sum += td.toll_factor * td.toll[i] + td.distance_factor * td.link_length[i]
    end
    return sum / norm_const
end

function mod_hessian(x::Vector{Float64}, td::NDFW.TA_Data)
    no_arc = length(td.init_node)

    h = zeros(no_arc, no_arc)
    h_diag = mod_hessian_diag(x, td, norm_const)

    for i in 1:no_arc
        h[i, i] = h_diag[i]
    end

    return h

    #Link travel time = free flow time * ( 1 + b * (flow/capacity)^Power ).
end

function mod_hessian_diag(x::Vector{Float64}, td::NDFW.TA_Data, norm_const=1)
    h_diag = Array{Float64}(undef, length(x))
    for i in 1:length(x)
        if td.power[i] >= 1.0
            h_diag[i] = td.free_flow_time[i] * td.b[i] * td.power[i] * (x[i]^(td.power[i] - 1)) / (td.capacity[i]^td.power[i])
        else
            h_diag[i] = 0 # Some cases, power is zero.
        end
    end
    # h_diag = free_flow_time .* b .* power .* (x.^(power-1)) ./ (capacity.^power)

    return h_diag / norm_const
    #Link travel time = free flow time * ( 1 + b * (flow/capacity)^Power ).
end


function modified_ta_frank_wolfe(td::NDFW.TA_Data; method=:bfw, max_iter_no=2000, step=:exact, log=:off, tol=1e-3, norm_const = 1)

    setup_time = time()

    if log==:on
        println("-------------------------------------")
        println("Network Name: $(td.network_name)")
        println("Method: $method")
        println("Line Search Step: $step")
        println("Maximum Interation Number: $max_iter_no")
        println("Tolerance for AEC: $tol")
        #println("Number of processors: ", nprocs())
    end

    # preparing a graph
    graph = create_graph(td.init_node, td.term_node)
    link_dic = sparse(td.init_node, td.term_node, collect(1:td.number_of_links))

    setup_time = time() - setup_time

    if log==:on
        println("Setup time = $setup_time seconds")
    end


    n_links = td.number_of_links

    iteration_time = time()


    # Finding a starting feasible solution
    travel_time = mod_BPR(zeros(n_links), td, norm_const)
    x0 = all_or_nothing(travel_time, td, graph, link_dic)

    # Initializing variables
    xk = x0
    tauk = 0.0
    yk_FW = x0
    sk_CFW = yk_FW
    Hk_diag = Vector{Float64}(undef, n_links)

    dk_FW = Vector{Float64}(undef, n_links)
    dk_bar = Vector{Float64}(undef, n_links)
    dk_CFW = Vector{Float64}(undef, n_links)
    dk = Vector{Float64}(undef, n_links)

    alphak = 0.0
    Nk = 0.0
    Dk = 0.0

    tauk = 0.0
    is_first_iteration = false
    is_second_iteration = false

    sk_BFW = yk_FW
    sk_BFW_old = yk_FW

    dk_bbar = Vector{Float64}(undef, n_links)
    muk = Vector{Float64}(undef, n_links)
    nuk = Vector{Float64}(undef, n_links)
    beta0 = 0.0
    beta1 = 0.0
    beta2 = 0.0

    trajectory = []

    for k in 1:max_iter_no
        # Finding yk

        travel_time = mod_BPR(xk, td, norm_const)
        yk_FW = all_or_nothing(travel_time, td, graph, link_dic) 

        trajectory = enter_trajectory_data(trajectory, xk, yk_FW, td, k, iteration_time, norm_const)

        # Basic Frank-Wolfe Direction
        dk_FW = yk_FW - xk 
        Hk_diag = mod_hessian_diag(xk, td, norm_const) # Hk_diag is a diagonal vector of matrix Hk

        # Finding a feasible direction
        if method == :fw # Original Frank-Wolfe
            dk = dk_FW 
        elseif method == :cfw # Conjugate Direction F-W
            if k==1 || tauk > 0.999999 # If tauk=1, then start the process all over again.
                sk_CFW = yk_FW 
                dk_CFW = sk_CFW - xk 
            else
                dk_bar = sk_CFW - xk  # sk_CFW from the previous iteration k-1
    
                Nk = dot( dk_bar, Hk_diag .* dk_FW ) 
                Dk = dot( dk_bar, Hk_diag .* (dk_FW - dk_bar) ) 
                
                delta = 0.0001
                if Dk != 0.0 && 0.0 <= Nk/Dk <= 1.0 - delta
                    alphak = Nk/Dk
                elseif Dk != 0.0 && Nk/Dk > 1.0 - delta
                    alphak = 1.0 - delta
                else
                    alphak = 0.0
                end
    
                # Generating new sk_CFW and dk_CFW
                sk_CFW = alphak .* sk_CFW .+ (1.0 - alphak) .* yk_FW
                dk_CFW = sk_CFW .- xk
            end
    
            # Feasible Direction to Use for CFW
            dk = dk_CFW
        elseif method == :bfw # Bi-Conjugate Direction F-W

            if tauk > 0.999999
                is_first_iteration = true
                is_second_iteration = true
            end

            if k==1 || is_first_iteration       # First Iteration is like FW
                sk_BFW_old = yk_FW
                dk_BFW = dk_FW
                is_first_iteration = false
            elseif k==2 || is_second_iteration  # Second Iteration is like CFW
                dk_bar = sk_BFW_old - xk # sk_BFW_old from the previous iteration 1

                Nk = dot( dk_bar, Hk_diag .* dk_FW )
                Dk = dot( dk_bar, Hk_diag .* (dk_FW - dk_bar) )

                delta = 0.0001
                if Dk != 0.0 && 0.0 <= Nk/Dk <= 1.0 - delta
                    alphak = Nk/Dk
                elseif Dk != 0.0 && Nk/Dk > 1.0 - delta
                    alphak = 1.0 - delta
                else
                    alphak = 0.0
                end

                # Generating new sk_BFW and dk_BFW
                sk_BFW = alphak .* sk_BFW_old .+ (1-alphak) .* yk_FW
                dk_BFW = sk_BFW .- xk

                is_second_iteration = false
            else
                # println("over there $tauk")
                # sk_BFW, tauk is from iteration k-1
                # sk_BFW_old is from iteration k-2

                dk_bar  = sk_BFW - xk
                dk_bbar = tauk * sk_BFW - xk + (1.0 - tauk) * sk_BFW_old

                muk = - dot( dk_bbar, Hk_diag .* dk_FW ) / dot( dk_bbar, Hk_diag .* (sk_BFW_old - sk_BFW) )
                nuk = - dot( dk_bar, Hk_diag .* dk_FW ) / dot( dk_bar, Hk_diag .* dk_bar) + muk*tauk/(1-tauk)

                muk = max(0.0, muk)
                nuk = max(0.0, nuk)

                beta0 = 1.0 / ( 1.0 + muk + nuk )
                beta1 = nuk * beta0
                beta2 = muk * beta0

                sk_BFW_new = beta0 * yk_FW + beta1 * sk_BFW + beta2 * sk_BFW_old
                dk_BFW = sk_BFW_new - xk

                sk_BFW_old = sk_BFW
                sk_BFW = sk_BFW_new
            end

            # Feasible Direction to Use for BFW
            dk = dk_BFW
        else
            error("The type of Frank-Wolfe method is specified incorrectly. Use :fw, :cfw, or :bfw.")
        end
        # dk is now identified.

        if step == :exact
            # Line Search from xk in the direction dk
            optk = optimize(tau -> mod_objective(xk + tau * dk, td, norm_const), 0.0, 1.0, GoldenSection())
            tauk = optk.minimizer
        elseif step == :newton
            # Newton step
            tauk = -dot(mod_gradient(xk, td, norm_const), dk) / dot(dk, Hk_diag .* dk)
            tauk = max(0.0, min(1.0, tauk))
        end

        # Average Excess Cost
        average_excess_cost = norm_const * ( xk' * travel_time - yk_FW' * travel_time ) / sum(td.travel_demand) 
        if log==:on
            @printf("k=%4d, tauk=%15.10f, objective=%15f, aec=%15.10f\n", k, tauk, mod_objective(xk, td, norm_const), average_excess_cost)
        end

        # rel_gap = ( objective(xk) - best_objective ) / best_objective

        # Convergence Test
        if average_excess_cost < tol
        # if rel_gap < tol
            break
        end

        # Update x
        new_x = xk + tauk * dk
        xk = new_x

        
        @assert minimum(xk) >= 0

    end


    iteration_time = time() - iteration_time

    if log==:on
        println("Iteration time = $iteration_time seconds")
    end

    return xk, travel_time, mod_objective(xk, td, norm_const), trajectory

end

function enter_trajectory_data(trajectory, xk, yk_FW, td, k, iteration_time, norm_const)
    g = mod_gradient(xk, td, norm_const)
    dual_gap = dot(g, xk - yk_FW)
    state = State(k, xk, mod_objective(xk, td, norm_const), dual_gap, time() - iteration_time)
    push!(trajectory, state)

    return trajectory
end
