using LinearAlgebra
using DifferentialEquations
using Distributions
using Random
using Printf
using Plots
using RCall
using CPUTime
using DynamicalSystems
using StatsPlots
# using StatsPlots
include("utils.jl")
include("lorenz_ode.jl")
include("abc.jl")
include("abcMCMC.jl")
include("BayesianCalibration.jl")
include("abcSMC.jl")


# sigma, rho, beta, x0
# sigma_prior = Uniform(5,15)
# rho_prior = Uniform(20, 30)
# beta_prior = Uniform(8/3 - 1.5, 8/3 +1.5)
# x0_prior = Uniform(7, 13)

# sigma_prior = truncated(Normal(10., 5), 0, 20)
# rho_prior = truncated(Normal(26, 10), 0, 40)
# beta_prior = truncated(Normal(8/3., 2), 1, 4)


sigma_prior = Uniform(5,15)
rho_prior = Uniform(20, 30)
beta_prior = Uniform(8/3 - 1.5, 8/3 +1.5)

sd = 0.01
sigma_true = truncated(Normal(10., sd), 0, Inf)
rho_true = truncated(Normal(26., sd), 0, Inf)
beta_true = truncated(Normal(8/3., sd), 0, Inf)

function data_generator_lorenz(p)

    u0 = [10., 10., 10.]
    time_window = (0., 20.)

    y = solve_lorenz(u0,time_window,p)
end

function data_generator_lorenz10(p)

    u0 = [10., 10., 10.]
    time_window = (0., 10.)

    y = solve_lorenz(u0,time_window,p)
end

true_prior=(sigma_true, rho_true, beta_true)
true_p = map(x -> rand(x),true_prior)
y = data_generator_lorenz(true_p)

y1 = data_generator_lorenz((10., 26., 8/3))
y2 = data_generator_lorenz((10, 26, 8/3 - 1e-3))

sol1=hcat(y1[:]...)'
sol2=hcat(y2[:]...)'
l2_error = sqrt(sum(sum((sol1.-sol2).^2,dims=2)*0.01,dims=1))[1]

#
#
#
# ds_1 = Systems.lorenz([10., 10., 10.];σ= 10., ρ= 26., β= 8/3)
# ds_2 = Systems.lorenz([10., 10., 10.];σ= 10. , ρ= 26., β= 8/3 + 1)
# lyap_1 = ChaosTools.lyapunov(ds_1, 1000.0, dt=1.0, Ttr=10.0)
# lyap_2 = ChaosTools.lyapunov(ds_2, 1000.0, dt=1.0, Ttr=10.0)
# dist_pert = abs(lyap_1 - lyap_2)

#abc
algo_parameter_ABC = (prior = ( sigma_prior, rho_prior, beta_prior),epsilon = 45.,
eta= identity_mapping, d= compute_l2_norm_decay)

#  ABC MCMC
# algo_parameter_mcmc = (prior = ( sigma_prior, rho_prior, beta_prior) ,epsilon = 60., eta = identity_mapping,
#                         d= compute_l2_norm_decay, proposal = proposal_Normal,
#                         proposalRatio = proposalRatio_Normal,sd = (0.5,0.5,0.5),
#                         thinning = 100, burn_in = 100,verbose=true)

algo_parameter_mcmc = (prior = ( sigma_prior, rho_prior, beta_prior) ,epsilon = 60., eta = identity_mapping,
                        d= compute_l2_norm_decay, proposal = proposal_Normal,
                        proposalRatio = proposalRatio_Normal,sd = (0.5,0.5,0.5),
                        thinning = 20, burn_in = 100,verbose=true)

#smc
algo_parameter_smc = (prior = ( sigma_prior, rho_prior, beta_prior), time_final=8, eps_list = [80,70, 65,60,55,50,44,45],
eta= identity_mapping, d= compute_l2_norm_decay, kernel=proposal_Normal,
kernel_density=proposal_Normal_density,
sd=(0.3,0.3, 0.3), resample_method=systematic_resample, verbose=true)

output_abc = ABC(y, data_generator_lorenz, algo_parameter_ABC, 250)
output_smc = ABC_SMC(y, data_generator_lorenz, algo_parameter_smc, 250)
output_mcmc = ABC_MCMC(y, data_generator_lorenz, algo_parameter_mcmc, 250)

n=100

sigma_abc=output_abc[1][:,1]; rho_abc=output_abc[1][:,2]; beta_abc=output_abc[1][:,3]

ds_observed = Systems.lorenz([10., 10., 10.];σ= true_p[1], ρ= true_p[2], β= true_p[3])
lyap_observed = ChaosTools.lyapunov(ds_observed, 1000.0, dt=10.0, Ttr=100.0)
dist_lyapunov_abc=zeros(n)
lya_abc=zeros(n)
for j=1:n
    u0 = [10., 10., 10.]
    ds_fitted = Systems.lorenz(u0;σ= sigma_abc[j], ρ= rho_abc[j], β= beta_abc[j])
    lyap_fit = ChaosTools.lyapunov(ds_fitted, 1000.0, dt=10.0, Ttr=100.0)
    dist_lyapunov_abc[j] = abs(lyap_fit-lyap_observed)/lyap_observed
    lya_abc[j] = lyap_fit
end

sigma_mcmc=output_mcmc[1][:,1]; rho_mcmc=output_mcmc[1][:,2]; beta_mcmc=output_mcmc[1][:,3]
dist_lyapunov_mcmc=zeros(n)
lya_mcmc=zeros(n)
for j=1:n
    u0 = [10., 10., 10.]
    ds_fitted = Systems.lorenz(u0;σ= sigma_mcmc[j], ρ= rho_mcmc[j], β= beta_mcmc[j])
    lyap_fit = ChaosTools.lyapunov(ds_fitted, 1000.0, dt=1.0, Ttr=10.0)
    dist_lyapunov_mcmc[j] = abs(lyap_observed-lyap_fit)/lyap_observed
    lya_mcmc[j]=lyap_fit
end

sigma_smc=output_smc[1][:,1]; rho_smc=output_smc[1][:,2]; beta_smc=output_smc[1][:,3]
dist_lyapunov_smc=zeros(n)
lyapunov_smc=zeros(n)
for j=1:n
    u0 = [10., 10., 10.]
    ds_fitted = Systems.lorenz(u0;σ= sigma_smc[j], ρ= rho_smc[j], β= beta_smc[j])
    lyap_fit = ChaosTools.lyapunov(ds_fitted, 1000.0, dt=10.0, Ttr=10.0)
    dist_lyapunov_smc[j] = abs(lyap_fit-lyap_observed)/lyap_observed
    lyapunov_smc[j] = lyap_fit
end
using LaTeXStrings
lyap_dist_plot = plot(plot(1:n,dist_lyapunov_abc, seriestype=:scatter, markersize=3, alpha=0.5, label="ABC"),
plot(1:n,dist_lyapunov_mcmc,seriestype=:scatter, markersize=3, alpha=0.5, label="MCMC"),
plot(1:n,dist_lyapunov_smc,seriestype=:scatter,  markersize=3, alpha=0.5,label="SMC"),
ylim=(0,1.1), xlabel=L"n", ylabel=L"$|\lambda^* - \lambda^n|$",
legend=:outertopright,
 dpi=350, layout=(3,1),
 foreground_color_legend=nothing)
plot(lya_abc,zeros(n),seriestype=:scatter)
plot!([lyap_observed],[0], seriestype=:scatter, color="red")

lyap_density_plot = density(lya_abc); vline!([lyap_observed])


n= size(output_abc[1])[1]
p_abc=plot(y, vars=(1,3),w=1,legend=false,title=L"\textrm{ABC}", dpi=300)
abc_best_ind = sortperm(dist_lyapunov_abc)
abc_best_ind = randperm(n)
# abc_best_ind = 1:1:n
output_abc_best = output_abc[1][abc_best_ind,:]
for k in 1:20
    p_pred = (output_abc_best[k,1], output_abc_best[k,2], output_abc_best[k,3])
    y_fitted = data_generator_lorenz(p_pred)
    plot!(p_abc,y_fitted, vars=(1,3),alpha=0.3, color="#BBBBBB", legend=false, xlim=(-20, 20),xlabel=L"x", ylabel=L"z")
end
p_abc


# mcmc_best_ind = sortperm(dist_lyapunov_mcmc)
mcmc_best_ind =randperm(n)
output_mcmc_best = output_mcmc[1][mcmc_best_ind,:]
p_mcmc=plot(y, vars=(1,3),w=1,legend=false,title=L"\textrm{ABC-MCMC}", dpi=300)
for k in 1:20
    p_pred = (output_mcmc_best[k,1], output_mcmc_best[k,2], output_mcmc_best[k,3])
    y_fitted = data_generator_lorenz(p_pred)
    plot!(p_mcmc,y_fitted, vars=(1,3),alpha=0.3, color="#BBBBBB", legend=false, xlim=(-20, 20),xlabel=L"x", ylabel=L"z")
end
p_mcmc

smc_best_ind = sortperm(dist_lyapunov_smc)
mcmc_best_ind =randperm(n)
# smc_best_ind = 1:1:n
output_smc_best = output_smc[1][smc_best_ind,:]
p_smc=plot(y, vars=(1,3),w=1,legend=false,title=L"\textrm{ABC-SMC}",dpi=300)
for k in 1:20
    p_pred = (output_smc_best[k,1], output_smc_best[k,2], output_smc_best[k,3])
    y_fitted = data_generator_lorenz(p_pred)
    plot!(p_smc,y_fitted, vars=(1,3),alpha=0.3, color="#BBBBBB", legend=false, xlim=(-20, 20),xlabel=L"x", ylabel=L"z")
end
p_smc

sigma_density=density(output_abc[1][:,1], xlabel="σ", label="ABC")
density!(output_mcmc[1][:,1], label="MCMC")
density!(output_smc[1][:,1], label="SMC")


rho_density=density(output_abc[1][:,2], xlabel="ρ", label="ABC")
density!(output_mcmc[1][:,2], label="MCMC")
density!(output_smc[1][:,2], label="SMC")


b_density=density(output_abc[1][:,3], xlabel="b", label="ABC")
density!(output_mcmc[1][:,3], label="MCMC")
density!(output_smc[1][:,3], label="SMC")
