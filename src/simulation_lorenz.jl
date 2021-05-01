using LinearAlgebra
using DifferentialEquations
using Distributions
using Random
using Printf
using Plots
using RCall
using CPUTime
using JLD
# using StatsPlots
include("utils.jl")
include("lorenz_ode.jl")
include("abc.jl")
include("abcMCMC.jl")
include("abcSMC.jl")
include("BayesianCalibration.jl")


# sigma, rho, beta, x0


#
# sigma_prior = Truncated(Normal(10., 5), 0, Inf)
# rho_prior = Truncated(Normal(26, 4), 0, Inf)
# beta_prior = Truncated(Normal(8/3., 1), 0, Inf)

# sigma_prior = truncated(Normal(10., 5), 0, 20)
# rho_prior = truncated(Normal(26, 10), 0, 40)
# beta_prior = truncated(Normal(8/3., 1), 1, 4)

sigma_prior = Truncated(Normal(10., 5), 0, Inf)
rho_prior = Truncated(Normal(26., 10), 0, Inf)
beta_prior = Truncated(Normal(8/3., 1), 0, Inf)

#
# sigma_prior = Uniform(5,15)
# rho_prior = Uniform(20, 30)
# beta_prior = Uniform(8/3 - 1.5, 8/3 +1.5)

sd = 0.5
# sigma_true = Truncated(Normal(10., sd), 1, 15)
# rho_true = Truncated(Normal(26., sd), 15, 40)
# beta_true = Truncated(Normal(8/3., sd_beta), 1, 4)
sigma_true = truncated(Normal(10., sd), 0, Inf)
rho_true = truncated(Normal(26., sd), 0, Inf)
beta_true = truncated(Normal(8/3., sd), 0, Inf)
x0_true = truncated(Normal(10., sd), 0, Inf)

# sigma_prior = Uniform(9,11)
# rho_prior = Uniform(25, 27)
# beta_prior = Uniform(8/3 - 0.5, 8/3 + 0.5)

true_prior=(sigma_true, rho_true, beta_true)


#abc
algo_parameter_ABC = (prior = ( sigma_prior, rho_prior, beta_prior),epsilon = 60.,
eta= identity_mapping, d= compute_l2_norm_decay)

#  ABC MCMC
algo_parameter_mcmc = (prior = ( sigma_prior, rho_prior, beta_prior) ,epsilon = 70., eta = identity_mapping,
                        d= compute_l2_norm_decay, proposal = proposal_Normal,
                        proposalRatio = proposalRatio_Normal,sd = (0.8,0.8,0.8),
                        thinning = 20, burn_in = 100,verbose=true)

# [ 65,62,60,58,55,52,48]
algo_parameter_smc = (prior = ( sigma_prior, rho_prior, beta_prior), time_final=5, eps_list = [ 80,70,65,62,60],
eta= identity_mapping, d= compute_l2_norm_decay, kernel=proposal_Normal,
kernel_density=proposal_Normal_density,
sd=(0.3,0.3, 0.3), resample_method=systematic_resample, verbose=true)


# Bayesian Calibration
algo_list = (ABC=ABC, ABC_MCMC=ABC_MCMC, ABC_SMC=ABC_SMC)

algo_param_list = (ABC = algo_parameter_ABC,
                   ABC_MCMC = algo_parameter_mcmc,
                   ABC_SMC = algo_parameter_smc)

@printf("Running %d simulations \n", 2^8)
simulation = BayesianCalibration(Int(2^8),Int(250),0.10,algo_list,algo_param_list,
ode_model = solve_lorenz,
initial_state = [10., 10., 10.],
time_window=(0,20.0),
add_noise = false,
true_p_dist=true_prior,
save_param_output=true)



# calibration
show_calibration(simulation)
# #
results_dir = "../results/simulation_lorenz5/"

save(results_dir*"output.jld","output",simulation)

using DataFrames
# simulation = load(results_dir*"output.jld")["output"]

function to_dataframe(x)
    tmp = DataFrame(x)
    rename!(tmp,[:ABC,:MCMC,:SMC])
end

CSV.write(results_dir*"calibration.csv",to_dataframe(simulation[2]))
CSV.write(results_dir*"cpu.csv",to_dataframe(simulation[4]))
CSV.write(results_dir*"ess.csv",to_dataframe(simulation[5]))
CSV.write(results_dir*"ess_cpu.csv",to_dataframe(simulation[6]))
