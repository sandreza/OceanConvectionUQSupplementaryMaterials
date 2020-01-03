include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")

use_covariance_estimate = true
case_range = 1:2
resolution_range = 1:1
const ensemble_size = 10^3
for resolution in resolutions[resolution_range]
    for case in cases[case_range]
        # construct filename
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        # construct default loss function
        N = resolution[1]
        Δt = resolution[2]
        ℒ = closure_default_loss_function(filename, N = N, Δt = Δt)
        # choose default parameters
        optimal_𝑪 = copy(default_𝑪)
        # optimize using default optimize, iterate a few times
        println("-------------------")
        println("For case $case ")
        println("and resolution " * string(resolution[1]))
        println("starting mcmc")
        Random.seed!(1234)
        resolution_label = "_res_" * string(resolution[1])
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_optima.jld2"
        mcmc_data = jldopen(filename, "r")
        initial_𝑪 = mcmc_data["parameter"]
        ℒ⁰ = mcmc_data["loss"]
        if use_covariance_estimate
            Σ = mcmc_data["covariance"]
        end
        close(mcmc_data)
        # scale the loss function by ℒ
        nll(𝑪) = ℒ(𝑪) / ℒ⁰
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
        # parameters for mcmc
        nt = ensemble_size
        frequency = 100
        # define proposal matrix, 5% of default value
        proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
        if use_covariance_estimate
            proposal = CoreFunctionality.closure_proposal(Σ, left_bounds = left_bounds, right_bounds = right_bounds)
        end
        # now markov chain
        CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename, verbose = true)
        println("done")
    end
end

###
# run mcmc on prior distribution, produces uniform distribution
nll(𝑪) = 1.0
filename = pwd() * "/mcmc_data/" * "prior"* "_mcmc.jld2"
# parameters for mcmc
nt = 20000
frequency = 100
# define proposal matrix, 5% of default value
σ = default_𝑪 * 0.5
σ[4] = 3
proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
# now markov chain
initial_𝑪 = copy(default_𝑪)
CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename)
println("done")
