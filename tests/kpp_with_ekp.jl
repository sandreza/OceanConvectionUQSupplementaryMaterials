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

nt = length(subsample)
coarseT = zeros(N, nt)
for i in 1:nt
    coarseT[:,i] .= CoreFunctionality.avg( les.T[:, subsample[i]], N)
end
# define the loss function
ℒ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=false, power = 2, f1 = mean, f2 = maximum )

# h 
params = [0.1, 6.33, 1.36, 3.19]
ℒ(params)

# G(params) = [ℒ(params)]
G(params) = 𝒢(params)[:]
##  EKP part
rng_seed = 41
Random.seed!(rng_seed)

# Number of synthetic observations from G(u)
n_obs = length(coarseT)
# Defining the observation noise level
noise_level =  1e-3   
# Independent noise for synthetic observations       
Γy = noise_level * Matrix(I, n_obs, n_obs) 
noise = MvNormal(zeros(n_obs), Γy)

# Loss Function Minimum
# y_obs  = Γy[1,:]
y_obs = coarseT[:]
# Define Prior
bounds = [(0,1), (0,10), (0,8), (0,6)]
prior_distns = [Parameterized(Uniform(bound...)) for bound in bounds]
constraints = [[bounded(bound...)] for bound in bounds]
# contraints = [[no_constraint()] for bound in bounds]
prior_names = ["cs", "nl", "kap", "ch"]
prior = ParameterDistribution(prior_distns, constraints, prior_names)
prior_mean = reshape(get_mean(prior),:)
prior_cov = get_cov(prior)

# Calibrate
N_ens = 100  # number of ensemble members
N_iter = 2 # number of EKI iterations
initial_ensemble = construct_initial_ensemble(prior, N_ens; rng_seed=rng_seed)
initial_losses = [G(initial_ensemble[i,:])[1] for i in 1:N_ens]
current_noise = minimum(initial_losses)
parameters = initial_ensemble
# y_obs .= current_noise
obs_mean = y_obs
obs_noise_cov = Γy
ekiobj = EnsembleKalmanProcess(initial_ensemble, y_obs, Γy, Inversion(), Δt = 1e-2)
##
g_ens = zeros(size(ekiobj.u[1])[1], size(G(ekiobj.u[1][1,:]))[1])

tic = time()
run_eki!(ekiobj, N_iter, G, g_ens, N_ens, prior)
toc = time()
println(toc-tic)

##
gr(size = (300,300))
index1 = 1
index2 = 2
for i in eachindex(ekiobj.u)
    p = plot(ekiobj.u[i][:,index1], ekiobj.u[i][:,index2], seriestype=:scatter, xlims = bounds[index1], ylims = bounds[index2])
    display(p)
    sleep(0.1)
end

##
params = ekiobj.u[end]
before_fail = ekiobj.u[12]
## start of fail 
check = ekiobj.u[12]
minimum(check[:])

##
scatter(initial_ensemble[:,1], log10.(initial_losses[:]), label = false)
ind1 = 1
ind2 = 2
scatter(initial_ensemble[:,ind1][:], initial_ensemble[:,ind2][:], log10.(initial_losses[:]))
