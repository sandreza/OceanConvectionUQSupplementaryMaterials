include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")


nonlocal_𝑪 = [0.09141046320860055, 0.0, 1.6322784658475666, 4.230574358293718];
#perform mcmc with no nonlocal term
for case in cases[1:1]
    for resolution in resolutions[1:1]
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
        initial_𝑪 = nonlocal_𝑪
        initial_𝑪[2] = 0.0
        ℒ⁰ = ℒ(initial_𝑪)
        # scale the loss function by ℒ
        nll(𝑪) = ℒ(𝑪) / ℒ⁰
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_no_nonlocal_mcmc.jld2"
        # parameters for mcmc
        nt = 20000
        frequency = 100
        # define proposal matrix, 5% of default value
        left_bounds = [0.0, 0.0, 0.0, 0.0]
        right_bounds = [1.0, eps(100.0), 6.0, 12.0]
        proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
        # now markov chain
        CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename)
        println("done")
    end
end
