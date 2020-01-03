include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")


use_covariance_estimate = true

#σ = default_𝑪 * 0.1
#σ[6] = 0.025
case_range = 3:1:34
case_range = [1,6,7,8,9,10]
case_range = 1:2
for resolution in resolutions[1:1]
    for case in cases[case_range]
        # construct filename
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        # construct default loss function
        les = CoreFunctionality.OceananigansData(filename)
        N = resolution[1]
        Δt = resolution[2]
        ℒ = closure_flexible_loss_function(filename, N = N, Δt = Δt, power = 1.0)
        # choose default parameters
        optimal_𝑪 = copy(default_𝑪)
        # optimize using default optimize, iterate a few times
        println("-------------------")
        println("For case $case ")
        println("and resolution " * string(resolution[1]))
        println("starting mcmc")
        Random.seed!(1234)
        extra_label = "_flexible_new"
        resolution_label = "_res_" * string(resolution[1])
        filename = pwd() * "/mcmc_data/" * case * resolution_label * extra_label *  "_optima.jld2"
        mcmc_data = jldopen(filename, "r")
        initial_𝑪 = mcmc_data["parameter"]
        ℒ⁰ = mcmc_data["loss"]
        if use_covariance_estimate
            Σ = mcmc_data["covariance"]
        end
        close(mcmc_data)
        # scale the loss function by ℒ
        nll(𝑪) = ℒ(𝑪) / ℒ⁰
        filename = pwd() * "/mcmc_data/" * case * resolution_label * extra_label *  "_mcmc.jld2"
        println(filename)
        # parameters for mcmc
        nt = 1000
        frequency = 100
        # define proposal matrix, 5% of default value
        NN = sqrt(les.α * les. g * les.bottom_T)
        σ = initial_𝑪 * 0.1
        #=
        left_bounds = [0.0, 3.0, 5.0, 0.0, 0.0, NN]
        right_bounds = [0.01, 5.0, 10.0, eps(1.0), 1.0, NN + eps(1.0)]
        =#
        left_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, NN]
        right_bounds = [0.01, 8.0, 10.0, eps(1.0), 1.0, NN + eps(1.0)]
        proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
        if use_covariance_estimate
            proposal = CoreFunctionality.closure_proposal(Σ, left_bounds = left_bounds, right_bounds = right_bounds)
        end
        # now markov chain
        CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename, verbose = true)
        println("done")
    end
end
