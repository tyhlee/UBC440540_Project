using LinearAlgebra
using DifferentialEquations
using Distributions
using Random
using Printf
using Plots
using RCall
using CPUTime
using StatsBase
using JLD
# using StatsPlots
include("utils.jl")
include("sir_ode.jl")
include("abc.jl")
include("abcMCMC.jl")
include("abcSMC.jl")
include("BayesianCalibration.jl")

# Bayesian Calibration
algo_list = (ABC=ABC, ABC_MCMC=ABC_MCMC, ABC_SMC=ABC_SMC)

true_prior = (truncated(Normal(2,0.5),0,Inf),truncated(Normal(1,0.25),0,Inf))

algo_parameter_ABC = (prior = true_prior,
                     epsilon = 0.1,eta = identity_mapping, d = compute_norm)

algo_parameter_mcmc = (prior = true_prior,
                        epsilon = 1, eta = identity_mapping,
                        d= compute_norm, proposal = proposal_Normal,
                        proposalRatio = proposalRatio_Normal,sd = (0.5,0.25),
                        thinning = 10, verbose=false,
                        burn_in = 100)


algo_parameters_smc = (prior = true_prior,
                        time_final=5, eps_list =[10, 5, 1, 0.5, 0.1],
                        eta= identity_mapping, d= distanceFunction,
                        kernel=proposal_Normal,
                        kernel_density=proposal_Normal_density,
                        sd=(0.5,0.25), resample_method=systematic_resample,
                        verbose=false)

algo_param_list = (ABC = algo_parameter_ABC,
                    ABC_MCMC = algo_parameter_mcmc,
                    ABC_SMC = algo_parameters_smc)

simulation = BayesianCalibration(Int(2^10),Int(250),0.10,algo_list,algo_param_list,
                                ode_model = solve_ode,
                                initial_state = [50000.0;4.5;350]./50000,
                                time_window=(0,20.0),
                                add_noise = false,
                                true_p_dist=true_prior)


y = solve_ode([50000.0;4.5;350]./50000,(0,20.0),rand.(true_prior))
plot(y[1,:])
plot!(y[2,:])
plot!(y[3,:])

# calibration
show_calibration(simulation)

results_dir = "results/simulation_SIR/"

save(results_dir*"output.jld","output",simulation)

using DataFrames
simulation = load(results_dir*"output.jld")["output"]

function to_dataframe(x)
    tmp = DataFrame(x)
    rename!(tmp,[:ABC,:MCMC,:SMC])
end

CSV.write(results_dir*"calibration.csv",to_dataframe(simulation[2]))
CSV.write(results_dir*"cpu.csv",to_dataframe(simulation[4]))
CSV.write(results_dir*"ess.csv",to_dataframe(simulation[5]))
CSV.write(results_dir*"ess_cpu.csv",to_dataframe(simulation[6]))
