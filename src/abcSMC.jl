# y: data
# yhat_generator: function that samples data given the parameters
# algo parameters: Namedtuple of algorithm parameters
# max_time: max time allowed for the algo
# N_samples: number of samples desired
function ABC_SMC(y,yhat_generator,algo_parameters,N_samples)
    CPUtic()
    resample_method = algo_parameters[:resample_method]
    kernel = algo_parameters[:kernel]
    kernel_density = algo_parameters[:kernel_density]
    sd = algo_parameters[:sd]
    eta_y = algo_parameters[:eta](y)
    size_y = size(y)
    q = length(algo_parameters[:prior])
    eps_list = algo_parameters[:eps_list]
    time_final = algo_parameters[:time_final]
    verbose = algo_parameters[:verbose]

    ess_list = zeros(time_final)
    ess_threshold = Int(N_samples*0.6)
    result = zeros(N_samples,q)
    weights = ones(N_samples)./N_samples
    theta = zeros(N_samples,q)
    t = 1
    # sample from the prior initially (same as rejection ABC)
    nParticle=1
    while nParticle < N_samples
        p = rand.(algo_parameters[:prior])
        yhat = yhat_generator(p)
        if algo_parameters[:d](eta_y,algo_parameters[:eta](yhat)) < eps_list[t]
            result[nParticle,:].=p
            nParticle +=1
        end
    end
    ess_list[t] = 1/(sum(weights.^2))
    t += 1

    result_prev = result
    weights_prev = weights
    while (t <= time_final)
        nParticle=1
        p_cand = Tuple(zeros(q))
        while nParticle <= N_samples
            # sample the index with with weights
            sampled_index = wsample(1:N_samples,weights_prev,1)
            p_star = Tuple(result[sampled_index,:])
            # propose candidate parameter from the kernel
            valid = true
            while valid
                p_cand = kernel(p_star, sd)
                valid = any(p_cand .< 0)
            end
            # if the density  of the prior at proposed param is 0, start again
            prior_cand = prod(pdf.(algo_parameters[:prior],p_cand))
            if prior_cand == 0
                continue
            end
            yhat = yhat_generator(p_cand)
            if (algo_parameters[:d](eta_y, algo_parameters[:eta](yhat)) < eps_list[t])
                # accept the proposed parameter
                result[nParticle,:] .= p_cand
                # Compute the weight of the current particle
                weight_denom = 0
                for j=1:N_samples
                    weight_denom += weights_prev[j]*prod((kernel_density.(Tuple(result_prev[j,:]), p_cand, sd) ))
                end
                weights[nParticle] = prior_cand/weight_denom
                nParticle+=1
            end

        end
        # Normalize weights
        weights = weights./sum(weights)
        # Compute effective sample size
        ess = 1/sum(weights.^2)
        ess_list[t] = ess
        # resample if the effective sample size is below the threshold
        if ess < ess_threshold
            idx = resample_method(weights)
            result = result[idx,:]
        end
        result_prev = result
        t += 1
    end
    t_end = CPUtoc_modified(false)
    if verbose
        @show ess_list
        @printf("Done SMC\n")
    end
    return result, t_end
end
