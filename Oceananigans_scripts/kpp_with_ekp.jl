include("../src/LocalOceanUQSupplementaryMaterials.jl")

include(pwd() * "/tests/ensemble_kalman_process.jl")

using Plots, Printf, Statistics, LinearAlgebra

# Get LES data
wd = pwd()
filename = wd * "/LES/general_strat_16_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

# Set parameters
N = 16
Δt = 10*60; # 10 minutes

# define the forward map
zᵖ = zeros(N)
start = argmin(abs.(les.t .- 86400 * 0.25)) #start at a quarter of a day
#start = 1
subsample = start:2:length(les.t)

# define the forward map
𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les,
                                              subsample = subsample, grid = zᵖ)




# define the loss function
ℒ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=false, power = 2, f1 = mean, f2 = maximum )

# h 
params = [0.1, 6.33, 3*1.36, 3.19*2]
ℒ(params)

##  EKP part

rng_seed = 41
Random.seed!(rng_seed)

# Number of synthetic observations from G(u)
n_obs = 1
# Defining the observation noise level
noise_level =  1e-8   
# Independent noise for synthetic observations       
Γy = noise_level * Matrix(I, n_obs, n_obs) 
noise = MvNormal(zeros(n_obs), Γy)

# Loss Function (unique minimum)
function G(u)
    return [sqrt((u[1]-1)^2 + (u[2]+1)^2)]
end

# Loss Function Minimum
u_star = [1.0, -1.0]
y_obs  = G(u_star) + 0 * rand(noise) 

# Define Prior
prior_distns = [Parameterized(Normal(0., sqrt(1))),
                Parameterized(Normal(-0., sqrt(1)))]
constraints = [[no_constraint()], [no_constraint()]]
prior_names = ["u1", "u2"]
prior = ParameterDistribution(prior_distns, constraints, prior_names)
prior_mean = reshape(get_mean(prior),:)
prior_cov = get_cov(prior)

# Calibrate
N_ens = 50  # number of ensemble members
N_iter = 20 # number of EKI iterations
initial_ensemble = construct_initial_ensemble(prior, N_ens;
                                                rng_seed=rng_seed)

ekiobj = EnsembleKalmanProcess(initial_ensemble,
                    y_obs, Γy, Inversion())
##
g_ens = zeros(size(ekiobj.u[1])[1], size(G(ekiobj.u[1][1,:]))[1])


tic = time()
run_eki!(ekiobj, N_iter, G, g_ens, N_ens)
toc = time()
println(tic-tic)