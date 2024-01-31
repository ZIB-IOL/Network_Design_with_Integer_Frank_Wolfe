function build_f_and_grad(alg, tn, gparams, cost_of_expansion, res, total_scenarios)
    # Impplement the objective function. This objective funciton is taken from the TransportationNetworks library. 

    od_pair_count = gparams.od_pair_count
    edge_count = length(gparams.init_nodes)
    new_edge_count = length(gparams.removed_edges)

    function f(x)
        if any(<(0), x) == true 
            if maximum(x[x .< 0]) < -1e-8 && (maximum(x[x .< 0]) > -1e-8)
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x[x .< 0] .= 0.0
        end

        total_scenario_vars = od_pair_count*edge_count + edge_count
        tot_sum = 0.0
        
        for scen in 1:total_scenarios
            first_var = (scen-1)*total_scenario_vars
            x_agg = @view(x[(first_var + od_pair_count*edge_count+1):(first_var + od_pair_count*edge_count+edge_count)])

            
            for i in eachindex(x_agg)
                tot_sum += tn.free_flow_time[i] * ( x_agg[i] + tn.b[i]* ( x_agg[i]^(tn.power[i]+1)) / (tn.capacity[i]^tn.power[i]) / (tn.power[i]+1))
                tot_sum += tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
                #sum += 0.0001*x_agg[i]
            end
        end

        tot_sum = tot_sum / total_scenarios

        

        x_nd = @view(x[(total_scenario_vars * total_scenarios +1):(total_scenario_vars * total_scenarios + new_edge_count)])
        for i in eachindex(x_nd)
            tot_sum += cost_of_expansion[i] * x_nd[i]
        end
        

        return tot_sum
    end

    #Gradient of objective function. 
    function grad!(storage, x)
        if any(x .< 0) == true 
            if maximum(x[x .< 0]) < -1e-8 && (maximum(x[x .< 0]) < -1e-8)
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x[x .< 0] .= 0.0
        end

        storage .= 0
        total_scenario_vars = od_pair_count*edge_count + edge_count
        for scen in 1:total_scenarios
            first_var = (scen-1)*total_scenario_vars
            x_agg = @view(x[(first_var + od_pair_count*edge_count+1):(first_var + od_pair_count*edge_count+edge_count)])
        
            @. storage[(first_var + od_pair_count*edge_count+1):(first_var + od_pair_count*edge_count+edge_count)] = (1.0/total_scenarios) * tn.free_flow_time * ( 1 + tn.b* ( x_agg^(tn.power)) / (tn.capacity^tn.power))# .+ 0.0001
        end

        x_nd = @view(x[(total_scenario_vars * total_scenarios+1):(total_scenario_vars * total_scenarios+new_edge_count)])


        
        @. storage[(total_scenario_vars * total_scenarios+1):(total_scenario_vars * total_scenarios + new_edge_count)] += cost_of_expansion
        

        #storage = (storage / res) 

        return storage

        #@info storage[(od_pair_count*edge_count+1):(od_pair_count*edge_count+edge_count)]

    end

    

    return f, grad!
end


function build_f_and_grad_indicator(alg, tn, gparams, cost_of_expansion, res, total_scenarios)
    # Impplement the objective function. This objective funciton is taken from the TransportationNetworks library. 

    od_pair_count = gparams.od_pair_count
    edge_count = length(gparams.init_nodes)
    new_edge_count = length(gparams.removed_edges)

    function f(x)
        if any(<(0), x) == true
            if maximum(x[x.<0]) < -1e-8 && (maximum(x[x.<0]) > -1e-8)
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x[x.<0] .= 0.0
        end

        total_scenario_vars = od_pair_count * edge_count + edge_count
        tot_sum = 0.0

        for scen in 1:total_scenarios
            first_var = (scen - 1) * total_scenario_vars
            x_agg = @view(x[(first_var+od_pair_count*edge_count+1):(first_var+od_pair_count*edge_count+edge_count)])


            for i in eachindex(x_agg)
                tot_sum += tn.free_flow_time[i] * (x_agg[i] + tn.b[i] * (x_agg[i]^(tn.power[i] + 1)) / (tn.capacity[i]^tn.power[i]) / (tn.power[i] + 1))
                tot_sum += tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
                #sum += 0.0001*x_agg[i]
            end
        end

        tot_sum = tot_sum / total_scenarios

        

        x_nd = @view(x[(total_scenario_vars*total_scenarios+1):(total_scenario_vars*total_scenarios+new_edge_count)])
        for i in eachindex(x_nd)
            tot_sum += cost_of_expansion[i] * (1 - x_nd[i])
        end
        

        return tot_sum
    end

    #Gradient of objective function. 
    function grad!(storage, x)
        if any(x .< 0) == true
            if maximum(x[x.<0]) < -1e-8 && (maximum(x[x.<0]) < -1e-8)
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x[x.<0] .= 0.0
        end

        storage .= 0
        total_scenario_vars = od_pair_count * edge_count + edge_count
        for scen in 1:total_scenarios
            first_var = (scen - 1) * total_scenario_vars
            x_agg = @view(x[(first_var+od_pair_count*edge_count+1):(first_var+od_pair_count*edge_count+edge_count)])

            @. storage[(first_var+od_pair_count*edge_count+1):(first_var+od_pair_count*edge_count+edge_count)] = (1.0 / total_scenarios) * tn.free_flow_time * (1 + tn.b * (x_agg^(tn.power)) / (tn.capacity^tn.power))# .+ 0.0001
        end

        x_nd = @view(x[(total_scenario_vars*total_scenarios+1):(total_scenario_vars*total_scenarios+new_edge_count)])


        
        @. storage[(total_scenario_vars*total_scenarios+1):(total_scenario_vars*total_scenarios+new_edge_count)] += -cost_of_expansion
        

        #storage = (storage / res)

        return storage

        #@info storage[(od_pair_count*edge_count+1):(od_pair_count*edge_count+edge_count)]

    end



    return f, grad!
end