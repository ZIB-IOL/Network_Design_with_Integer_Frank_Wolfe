function build_f_and_grad_penalty_method(alg, tn, gparams, cost_of_expansion, res, total_scenarios, use_bigm, penalty, pv=1.0)
    # Impplement the objective function. This objective funciton is taken from the TransportationNetworks library. 

    @assert pv > 1

    od_pair_count = gparams.od_pair_count
    edge_count = length(gparams.init_nodes)
    new_edge_count = length(gparams.removed_edges)
    removed_edge_idx = [gparams.link_dic[removed_edge[1], removed_edge[2]] for removed_edge in gparams.removed_edges]

    td_2_dest = [maximum([sum(gparams.travel_demand[s][:, z]) for s in 1:total_scenarios]) for z in 1:od_pair_count]

    function f(x)
        if any(<(0), x) == true
            if maximum(x[x.<0]) < -1e-8 && (maximum(x[x.<0]) > -1e-8)
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x[x.<0] .= 0.0
        end

        total_scenario_vars = od_pair_count*edge_count + edge_count
        totsum = 0.0
        
        # for scen in 1:total_scenarios
        #     first_var = (scen-1)*total_scenario_vars
        #     x_agg = @view(x[(first_var + od_pair_count*edge_count+1):(first_var + od_pair_count*edge_count+edge_count)])

        #     for i in eachindex(x_agg)
        #         totsum += tn.free_flow_time[i] * (x_agg[i] + tn.b[i] * (x_agg[i]^(tn.power[i] + 1)) / (tn.capacity[i]^tn.power[i]) / (tn.power[i] + 1))
        #         totsum += tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
        #         #totsum += 0.01 * x_agg[i]
        #     end
        # end

        for scen in 1:total_scenarios
            first_var = (scen - 1) * total_scenario_vars
            for i in 1:edge_count

                x_agg = sum(x[(first_var+(z-1)*edge_count+i)] for z in 1:od_pair_count)

                totsum += tn.free_flow_time[i] * (x_agg + tn.b[i] * (x_agg^(tn.power[i] + 1)) / (tn.capacity[i]^tn.power[i]) / (tn.power[i] + 1))
                totsum += tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
                #totsum += 0.01 * x_agg[i]
            end
        end

        totsum = totsum / total_scenarios

        if (alg == "IFW" || alg == "NLMO-IFW" || alg == "IFW-P" || alg == "NLMO-P" || alg == "BNDLMO-P")
            x_nd = @view(x[(total_scenario_vars * total_scenarios+1):(total_scenario_vars * total_scenarios+new_edge_count)])
            for i in eachindex(x_nd)
                totsum += cost_of_expansion[i] * x_nd[i]
            end
        end

        bigM = tn.total_od_flow

        # for scen in 1:total_scenarios
        #     first_var = (scen-1)*total_scenario_vars
        #     x_agg = @view(x[(first_var + od_pair_count*edge_count+1):(first_var + od_pair_count*edge_count+edge_count)])
        #     totsum += penalty * sum(max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0)^pv for j in eachindex(x_nd))
        # end

        @inbounds for scen in 1:total_scenarios
            first_var = (scen - 1) * total_scenario_vars
            #curr_x = @view(x[(first_var+1):(first_var+od_pair_count*edge_count)])

            for j in 1:new_edge_count
                if use_bigm == false
                    for z in 1:od_pair_count
                        curr_x = @view(x[(first_var+(z-1)*edge_count+1):(first_var+z*edge_count)])                    
                        totsum += penalty * max(curr_x[removed_edge_idx[j]] - td_2_dest[z] * x_nd[j], 0)^pv
                    end
                else
                    x_agg = sum(x[(first_var+(z-1)*edge_count+removed_edge_idx[j])] for z in 1:od_pair_count)
                    totsum += penalty * max(x_agg - bigM * x_nd[j], 0)^pv
                end
            end
        end

        return totsum
    end

    #Gradient of objective function. 
    function grad!(storage, x)
        if any(<(0), x)
            if maximum(x[x.<0]) < -1e-8 && (maximum(x[x.<0]) < -1e-8)
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x[x.<0] .= 0.0
        end

        storage .= 0

        total_scenario_vars = od_pair_count*edge_count + edge_count

        prob = (1.0/total_scenarios)

        x_nd = @view(x[(total_scenario_vars * total_scenarios+1):(total_scenario_vars * total_scenarios+new_edge_count)])

        @inbounds for scen in 1:total_scenarios
            first_var = (scen-1)*total_scenario_vars

            for i in 1:edge_count
                x_agg = sum(x[(first_var+(lz-1)*edge_count+i)] for lz in 1:od_pair_count)

                rhs = prob * tn.free_flow_time[i] * (1.0 + tn.b[i] * (x_agg^(tn.power[i])) / (tn.capacity[i]^tn.power[i]))
                
                for z in 1:od_pair_count

                    #x_agg = @view(x[(first_var + od_pair_count*edge_count+1):(first_var + od_pair_count*edge_count+edge_count)])
                    #curr_x = @view(x[first_var + (z-1)*edge_count + i])

                    storage[(first_var+(z-1)*edge_count+i)] += rhs 
                end
            end


            # g = zeros(edge_count)
            # h = zeros(new_edge_count)
            # bigM = tn.total_od_flow

            # for j in eachindex(x_nd)
            #     storage[first_var + od_pair_count*edge_count + removed_edge_idx[j]] += pv * penalty * (max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0)^(pv - 1))
            #     storage[total_scenario_vars * total_scenarios + j] += pv * penalty * (max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0)^(pv - 1)) * (-bigM)
            # end
        end

        if (alg == "IFW" || alg == "NLMO-IFW" || alg == "IFW-P" || alg == "NLMO-P" || alg == "BNDLMO-P")
            storage[(total_scenario_vars*total_scenarios+1):(total_scenario_vars*total_scenarios+new_edge_count)] += cost_of_expansion
        end

        bigM = tn.total_od_flow

        # for scen in 1:total_scenarios
        #     first_var = (scen - 1) * total_scenario_vars
        #     x_agg = @view(x[(first_var+od_pair_count*edge_count+1):(first_var+od_pair_count*edge_count+edge_count)])
        #     for j in eachindex(x_nd)
        #         storage[first_var+od_pair_count*edge_count+removed_edge_idx[j]] += pv * penalty * (max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0)^(pv - 1))
        #         storage[total_scenario_vars*total_scenarios+j] += pv * penalty * (max(x_agg[removed_edge_idx[j]] - bigM * x_nd[j], 0)^(pv - 1)) * (-bigM)
        #     end
        # end

        for scen in 1:total_scenarios
            first_var = (scen - 1) * total_scenario_vars
            for j in 1:new_edge_count
                if use_bigm == false
                    for z in 1:od_pair_count
                        curr_x = @view(x[(first_var+(z-1)*edge_count+1):(first_var+z*edge_count)])
                        storage[first_var+(z-1)*edge_count+removed_edge_idx[j]] += pv * penalty * (max(curr_x[removed_edge_idx[j]] - td_2_dest[z] * x_nd[j], 0)^(pv - 1))
                        storage[total_scenario_vars*total_scenarios+j] += pv * penalty * (max(curr_x[removed_edge_idx[j]] - td_2_dest[z] * x_nd[j], 0)^(pv - 1)) * (-td_2_dest[z])
                    end
                else
                    x_agg = sum([x[(first_var+(z-1)*edge_count+removed_edge_idx[j])] for z in 1:od_pair_count])
                    for z in 1:od_pair_count
                        storage[first_var+(z-1)*edge_count+removed_edge_idx[j]] += pv * penalty * (max(x_agg - bigM * x_nd[j], 0)^(pv - 1))
                    end
                    storage[total_scenario_vars*total_scenarios+j] += pv * penalty * (max(x_agg - bigM * x_nd[j], 0)^(pv - 1)) * (-bigM)
                end
            end
        end


        # storage[(od_pair_count*edge_count+1):(od_pair_count*edge_count+edge_count)] += g
        # storage[((od_pair_count+1)*edge_count+1):((od_pair_count+1)*edge_count+new_edge_count)] += h
        
        return storage

    end

    return f, grad!
end