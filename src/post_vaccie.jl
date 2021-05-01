using JLD
using DataFrames
using Dates
using LinearAlgebra
using DifferentialEquations
using Distributions
using Random
using CSV
using Printf
using RCall
using CPUTime
using StatsPlots
include("utils.jl")
include("sir_ode.jl")
include("abc.jl")
include("abcMCMC.jl")
include("abcSMC.jl")

df_post = DataFrame(CSV.File("../data/covid19_bc_post_vaccine.csv"))


 """ Post-vaccine analysis """
 tot_days = size(df_post)[1]
 time_interval = 7
 
 y_post = hcat([df_post[:S],df_post[:I],df_post[:R]]...)[1:7:tot_days,:]
 y_post = y_post./1e4
 y_dim = size(y_post)[1]

 function data_generator_week(p)
     initial_state = y_post[1,:]; time_window=(1,Float64(tot_days)/time_interval)
     solve_ode(initial_state,time_window,p)
 end
 
 algo_parameters_abc = (prior = (Truncated(Normal(0.1,0.5),0,Inf),Truncated(Normal(0.2,0.2),0,Inf)),epsilon = 0.1,
 eta= identity_mapping, d= distanceFunction)
 output_abc=ABC(y_post', data_generator_week, algo_parameters_abc, 500)
 
 algo_parameter_mcmc = (prior = (Truncated(Normal(0.1,0.5),0,Inf),Truncated(Normal(0.2,0.2),0,Inf)),epsilon = 1, eta = identity_mapping,
                         d= compute_norm, proposal = proposal_Normal,
                         proposalRatio = proposalRatio_Normal,sd = (0.25,0.25),
                         thinning = 10000, burn_in=100, verbose=true)
 
 output_mcmc = ABC_MCMC(y_post', data_generator_week, algo_parameter_mcmc, 500)
 
 
 algo_parameters_smc = (prior = (Truncated(Normal(0.1,0.5),0,Inf),Truncated(Normal(0.2,0.2),0,Inf)), time_final=5, eps_list = [10, 5, 1, 0.5, 0.1],
 eta= identity_mapping, d= distanceFunction, kernel=proposal_Normal, 
 kernel_density=proposal_Normal_density,
 sd=(0.5,0.5), resample_method=systematic_resample, verbose=true)
 output_smc = ABC_SMC(y_post', data_generator_week, algo_parameters_smc, 500)
 
 
 abc_solutions=zeros(3,y_dim,500)
 abc_mcmc_solutions=zeros(3,y_dim,500)
 abc_smc_solutions=zeros(3,y_dim,500)
 for i=1:500
     abc_solutions[:,:,i] = data_generator_week(output_abc[1][i,:])[1:3,:]
     abc_mcmc_solutions[:,:,i] = data_generator_week(output_mcmc[1][i,:])[1:3,:]
     abc_smc_solutions[:,:,i] = data_generator_week(output_smc[1][i,:])[1:3,:]
 end
 
 save("../results/abc_solution_500_post2.jld","abc_solutions",abc_solutions)
 save("../results/abc_mcmc_solution_500_post2.jld","abc_mcmc_solutions",abc_mcmc_solutions)
 save("../results/abc_smc_solution_500_post2.jld","abc_smc_solutions",abc_smc_solutions)
 
 alpha_level=1
 beta_plot = density(output_abc[1][:,1],alpha=alpha_level, label="ABC",  title="Predicted β Posterior", xlabel="β", ylabel="Density")
 density!(output_mcmc[1][:,1],alpha=alpha_level, label="ABC-MCMC")
 density!(output_smc[1][:,1],alpha=alpha_level, label="ABC-SMC")
 savefig(beta_plot,"../figs/beta_plot_post2.png")
 
 gamma_plot = density(output_abc[1][:,2],alpha=alpha_level, label="ABC",  title="Predicted γ Posterior", xlabel="γ", ylabel="Density")
 density!(output_mcmc[1][:,2],alpha=alpha_level, label="ABC-MCMC")
 density!(output_smc[1][:,2],alpha=alpha_level, label="ABC-SMC")
 savefig(gamma_plot,"../figs/gamma_plot_post2.png")
 
 # ess
ess_output = zeros(3)
ess_output_time = copy(ess_output)
ess_output[1]  = try rcopy(R"mean(mcmcse::ess($(output_abc[1])))") catch; 0 end
ess_output_time[1]  = ess_output[1] / output_abc[2]
ess_output[2]  = try rcopy(R"mean(mcmcse::ess($(output_mcmc[1])))") catch; 0 end
ess_output_time[2]  = ess_output[2] / output_mcmc[2]
ess_output[3]  = try rcopy(R"mean(mcmcse::ess($(output_smc[1])))") catch; 0 end
ess_output_time[3]  = ess_output[3] / output_smc[2]

save("../results/covid/ess_post_2.jld","ess_output",ess_output)
save("../results/covid/ess_cputime_post_2.jld","ess_output_time",ess_output_time)


