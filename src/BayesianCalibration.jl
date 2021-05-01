# outputs I want
# true parameter value list
#
# TODO: how to compute the posterior distribution of parameters in ABC
function BayesianCalibration(N_experiments,N_samples,alpha,
    algo_list,algo_parameter_list;
    ode_model = solve_ode,
    initial_state = [99.0;1.0;0.0],
    time_window=(0,10.0),
    add_noise = false,
    true_p_dist=(Gamma(2,1),Gamma(1,1)),
    save_param_output=false)

    # initial pop
    s0 = sum(initial_state)

    # num of algorithms
    N_algo = length(algo_list)

    # N param
    N_param = length(true_p_dist)

    # list of param values
    true_param = rand.(true_p_dist,N_experiments)

    # Calibration value
    res = zeros(N_experiments,N_algo)

    # algo output
    if save_param_output
        param_output = zeros(N_experiments,N_samples,N_param,N_algo)
    end
    # algo time
    cputime_output = zeros(N_experiments, N_algo)

    # ess
    ess_output = zeros(N_experiments,N_algo)

    # ess per unit time
    ess_output_time = copy(ess_output)

    # data generator
    function data_generator(p)
        y = ode_model(initial_state,time_window,p)
        if add_noise
            y = y + rand(LogNormal(0,0.5),size(y))
            y = y ./ sum(y,dims=1) * s0
        end
        y
    end


    for i in 1:N_experiments
        # generate data
        # firt sample true SIR parameters
        # true_p = (true_param[1][i],true_param[2][i])
        true_p_vec=zeros(N_param)
        for k in 1:N_param
            true_p_vec[k] = true_param[k][i]
        end
        true_p = Tuple(true_p_vec)
        # generate data
        y = data_generator(true_p)

        # iterate over algos
        for algo_index in 1:N_algo
            # perform algo
            output, cpu_time = algo_list[algo_index](y,data_generator,algo_parameter_list[algo_index],N_samples)
            # check whether parameters are contained in the 95% posterior dist
            # given independence both have to be contained

            res[i,algo_index] = all(map(ind -> check_in_interval(true_p[ind],compute_posterior(vec(output[:,ind]),alpha)),1:length(true_p)))

            # save output
            if save_param_output
                param_output[i,:,:,algo_index] = output
            end
            cputime_output[i,algo_index] =  cpu_time
            # use R's package mcmcse to compute ESS
            ess_output[i,algo_index]  = try rcopy(R"mean(mcmcse::ess($output))") catch; 0 end
            ess_output_time[i,algo_index]  = ess_output[i,algo_index] / cpu_time
        end

        if !save_param_output
            param_output = nothing
        end

        if i % 10 == 0
            @show(i)
        end
    end

    return(true_param = true_param,
           calibration = res,
           simulated_param = param_output,
           cpu_time = cputime_output,
           ess = ess_output,
           ess_time = ess_output_time,
           algo_param_list = algo_parameter_list)
end
