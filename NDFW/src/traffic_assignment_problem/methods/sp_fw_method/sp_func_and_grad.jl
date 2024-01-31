function build_f_and_grad_sp(alg, tn, gparams, norm_const)
    # Impplement the objective function. This objective funciton is taken from the TransportationNetworks library. 

    od_pair_count = gparams.od_pair_count
    edge_count = length(gparams.init_nodes)

    function f(x)
        if any(<(0), x) == true 
            if maximum(x[x .< 0]) < -1e-8 && (maximum(x[x .< 0]) > -1e-8)
                @info "Observed negative x with most negative being $(maximum(x[x.<0]))"
            end
            x[x .< 0] .= 0.0
        end

        tot_sum = 0.0
        
        
        for i in eachindex(x)
            tot_sum += tn.free_flow_time[i] * ( x[i] + tn.b[i]* ( x[i]^(tn.power[i]+1)) / (tn.capacity[i]^tn.power[i]) / (tn.power[i]+1))
            tot_sum += tn.toll_factor * tn.toll[i] + tn.distance_factor * tn.link_length[i]
        end

        return tot_sum / norm_const
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

    
        for i in 1:edge_count
            storage[i] = tn.free_flow_time[i] * (1 + tn.b[i] * (x[i]^(tn.power[i])) / (tn.capacity[i]^tn.power[i])) / norm_const
        end
        
        return storage


    end

    

    return f, grad!
end


