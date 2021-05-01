# y: data
# yhat_generator: function that samples data given the parameters
# algo parameters: Namedtuple of algorithm parameters
# N_samples: number of samples desired
function ABC_MCMC(y,yhat_generator,algo_parameters,N_samples)
    CPUtic()

    proposal = algo_parameters[:proposal]
    proposalRatio = algo_parameters[:proposalRatio]
    sd = algo_parameters[:sd]
    eta_y = algo_parameters[:eta](y)
    size_y = size(y)
    q = length(algo_parameters[:prior])
    thinning_interval = algo_parameters[:thinning]
    N_samples = N_samples + algo_parameters[:burn_in]
    result = zeros(N_samples,q)
    dist = Inf
    p_prev=Tuple(zeros(q))

    # initial point using ABC
    while (dist > algo_parameters[:epsilon])
        p_prev = rand.(algo_parameters[:prior])
        yhat = yhat_generator(p_prev)
        if size_y != size(yhat)
            continue
        end
        dist = algo_parameters[:d](eta_y, algo_parameters[:eta](yhat))
    end

    naccept=0
    k = 1
    i = 0
    p_cand = Tuple(zeros(q))
    while i <= (N_samples*thinning_interval)
        # propose candidate parameters in log space
        valid = true
        while valid
            p_cand = proposal(p_prev, sd)
            valid = any(p_cand .< 0)
        end
        # p_cand = proposal(p_prev, sd)


        # p_cand = kernel(p_prev, sd)[1]
        #generate data
        yhat = yhat_generator(p_cand)
        # yhat = yhat_generator( p_cand)
        if size_y != size(yhat)
            continue
        end
        i += 1

        dist = algo_parameters[:d](eta_y, algo_parameters[:eta](yhat))

        u = rand(Uniform(0,1))

        proposal_ratio = proposalRatio_Normal(p_prev,p_cand,sd)

        log_prior_prev = sum(logpdf.(algo_parameters[:prior], p_prev))
        log_prior_cand = sum(logpdf.(algo_parameters[:prior], p_cand))

        log_alpha = log_prior_cand + sum(log.(proposal_ratio)) - log_prior_prev
        if (i % thinning_interval==0)
            if (log(u) < log_alpha && dist < algo_parameters[:epsilon])
                result[k,:] .= p_cand
                p_prev = p_cand
                naccept += 1
            else
                result[k,:] .= p_prev
            end
            k += 1
        end
    end
    if algo_parameters[:verbose]
        @printf("acceptance rate=%f\n", naccept/N_samples)
    end
    if algo_parameters[:burn_in] > 0
        result = result[(1+algo_parameters[:burn_in]):end,:]
    end
    t_end = CPUtoc_modified(false)

    return result,t_end
end
