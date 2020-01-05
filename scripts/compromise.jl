include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")

# boolean label for optimizing or performing mcmc
optimize_compromise = true
mcmc_compromise = true
use_covariance_estimate = true
const ensemble_size = 10^3
case_name = "compromise"
# default resolution
resolution = resolutions[1]
N = resolution[1]
Δt = resolution[2]
resolution_label = "_res_" * string(resolution[1])

# first define loss functions
case = cases[1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
# construct default loss function
ℒ1 = closure_default_loss_function(filename, N = N, Δt = Δt)
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain1 = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
indmin1 = argmin(e1)
𝑪1 = chain1[:,indmin1]
close(mcmc_data)


# now second case
case = cases[2]
#
filename = pwd() * "/LES/" * case * "_profiles.jld2"
# construct default loss function
ℒ2 = closure_default_loss_function(filename, N = N, Δt = Δt)
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain2 = mcmc_data["𝑪"]
e2 = mcmc_data["ε"]
indmin2 = argmin(e2)
𝑪2 = chain2[:,indmin2]
close(mcmc_data)
# define loss function
optimal_𝑪 = ( 𝑪1 .+ 𝑪2 ) ./ 2
a = ℒ1(𝑪1)
b = ℒ1(𝑪2)
c = ℒ2(𝑪2)
d = ℒ2(𝑪1)

# now define combined loss function
scale = (c+d) / (a+b)
ℒ_compromise(𝑪) = 0.5 *( ℒ1(𝑪) * scale + ℒ2(𝑪) )

if optimize_compromise == true
    # get optimal parameters in order to properly scale contribution of loss functions
    default_𝑪 = ( 𝑪1 .+ 𝑪2 ) ./ 2
    println("starting optimization")
    Random.seed!(1234)
    default_ℒ = ℒ_compromise(default_𝑪)
    # random walk optimization
    println("random walk")
    # optimal_𝑪 = CoreFunctionality.optimize(default_𝑪, nll; nt = 1000, restart = 2, proposal = proposal, rescale = true, freq = 100, scale = default_ℒ)
    optimal_𝑪, Σ = CoreFunctionality.optimize_and_estimate_proposal(default_𝑪, ℒ_compromise, left_bounds, right_bounds, nt = 1000, restart = 2, proposal = [], filename = [], rescale = true, freq = 100, verbose = true)
    default_ℒ = ℒ_compromise(default_𝑪)
    optimal_ℒ = ℒ_compromise(optimal_𝑪)
    println("--------------------")
    println("The default parameters are $default_𝑪")
    println("The lossfunction value is $(default_ℒ)")
    println("The optimal parameters are $optimal_𝑪")
    println("The lossfunction value is $optimal_ℒ")
    println("The improvement is a factor of $(default_ℒ/optimal_ℒ)")
    println("The covariance estimate is ")
    display(Σ)
    println("-------------------")
    # save optimal values and loss function value
    resolution_label = "_res_" * string(resolution[1])
    filename = pwd() * "/mcmc_data/" * case_name * resolution_label * "_optima.jld2"
    parameter = optimal_𝑪
    loss = optimal_ℒ
    covariance = Σ
    @save filename parameter loss covariance
end

if mcmc_compromise == true
    filename = pwd() * "/mcmc_data/" * case_name * resolution_label * "_optima.jld2"
    mcmc_data = jldopen(filename, "r")
    initial_𝑪 = mcmc_data["parameter"]
    ℒ⁰ = mcmc_data["loss"]
    if use_covariance_estimate
        Σ = mcmc_data["covariance"]
    end
    close(mcmc_data)
    # scale the loss function by ℒ
    nll(𝑪) = ℒ_compromise(𝑪) / ℒ⁰
    filename = pwd() * "/mcmc_data/" * case_name * resolution_label * "_mcmc.jld2"
    # parameters for mcmc
    nt = ensemble_size
    frequency = 100
    # define proposal matrix, 5% of default value
    proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
    if use_covariance_estimate
        proposal = CoreFunctionality.closure_proposal(Σ, left_bounds = left_bounds, right_bounds = right_bounds)
    end
    # now markov chain
    CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename)
    println("done")
end
