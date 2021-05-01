function check_in_interval(num,interval)
    (num >= interval[1]) & (num <= interval[2])
end

function compute_posterior(vec_num,alpha)
    quantile(vec_num,[alpha/2,1-alpha/2])
end

function distanceFunction(ytilde, y)
    return mean((ytilde[3,:]-y[3,:]) .^ 2)
end

function identity_mapping(y)
    y
end

function compute_norm(y,yhat;p=2)
    norm(y[2:3,2:end] - yhat[2:3,2:end], p)/size(y)[1]
end

function compute_l2_norm(y,yhat)
    sol1=hcat(y[:]...)'
    sol2=hcat(yhat[:]...)'
    l2_error = sqrt(sum(sum((sol1.-sol2).^2,dims=2)*0.01,dims=1))[1]
end

function compute_l2_norm_decay(y,yhat)
    sol1=hcat(y[:]...)'
    sol2=hcat(yhat[:]...)'
    T=size(sol1)[1]
    l2_error = sqrt(sum(sum((sol1.-sol2).^2,dims=2).*exp.(-0.001*(1:1:T))*0.01,dims=1))[1]
end

function compute_full_norm(y,yhat;p=Inf)
    norm(y - yhat, p)/size(y)[1]
end

function compute_at_steps(y,yhat;index=3:3:18,p=Inf)
    norm((y .- yhat)[2:3,index], p)/length(index)
end



function random_walk(p_prev, sd)
    # return proposed p, log(q(p_n+1|p_n))

    p_propose = map( (x,y) -> x + rand((Normal(0,y))) ,p_prev,sd)

    log_pdf = map( (x,y,z) -> logpdf((Normal(x,z)), y ), p_prev, p_propose, sd )
    return p_propose, log_pdf
end


function proposal_Normal_density(p1, p2, sd)
    return map( (x,y,z) -> pdf(Normal(x,z),y), p1, p2, sd )
end

function proposal_LN(log_p_prev,sd)
    log_p_prev .+ rand.(Normal.(0,sd))
end

function proposalRatio_LN(p_prev,p_cand,sd)
    0.5
end

function proposal_Gamma(log_p_prev,k)
    #TODO
end

function proposalRatio_Gamma(p_prev,p_cand,k)
    # TODO
end

function proposal_Normal(p_prev,sd)
    p_prev .+ rand.(Normal.(0,sd))
end

function proposalRatio_Normal(p_prev,p_cand,sd)
    1
end

function systematic_resample(weights)
    N = length(weights)
    partitions = rand()/N .+ (0:1:N-1)./N
    idx = zeros(Int,N)
    Q = cumsum(weights)
    i=1;j=1;
    while i<=N
        if partitions[i] < Q[j]
            idx[i]=j
            i +=1
        else
            j+=1
        end
    end
    return idx
end

function CPUtoc_modified(verbose)
    t = CPUTime.CPUtoq()
    if verbose
        println("elpased CPU time: ",t," seconds")
    end
    return t
end

function show_calibration(simulation)
    @printf("Calibration: ABC %f, ABC MCMC %f, ABC SMC %f\n", mean(simulation[2],dims=1)[1], mean(simulation[2],dims=1)[2],mean(simulation[2],dims=1)[3])
    @printf("CPU_time: ABC %f, ABC MCMC %f, ABC SMC %f\n", mean(simulation[4],dims=1)[1], mean(simulation[4],dims=1)[2],mean(simulation[4],dims=1)[3])
    @printf("CPU_time: ABC %f, ABC MCMC %f, ABC SMC %f\n", median(simulation[4],dims=1)[1], median(simulation[4],dims=1)[2],median(simulation[4],dims=1)[3])
    @printf("ESS: ABC %f, ABC MCMC %f, ABC SMC %f\n", median(simulation[5],dims=1)[1], mean(simulation[5],dims=1)[2],mean(simulation[5],dims=1)[3])
    @printf("ESS: ABC %f, ABC MCMC %f, ABC SMC %f\n", mean(simulation[5],dims=1)[1], median(simulation[5],dims=1)[2],median(simulation[5],dims=1)[3])
    @printf("ESS Time: ABC %f, ABC MCMC %f, ABC SMC %f\n", mean(simulation[6],dims=1)[1], mean(simulation[6],dims=1)[2],mean(simulation[6],dims=1)[3])
    @printf("ESS Time: ABC %f, ABC MCMC %f, ABC SMC %f\n", median(simulation[6],dims=1)[1], median(simulation[6],dims=1)[2],median(simulation[6],dims=1)[3])
end
