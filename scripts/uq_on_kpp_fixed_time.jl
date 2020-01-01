include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")


case_range = 1:1
resolution_range = 1:1
days = [1.0, 2.0 , 3.0 , 4.0 , 5.0, 6.0, 7.0, 8.0]

for resolution in resolutions[resolution_range]
    for case in cases[case_range]
        for day in days
            # construct filename
            filename = pwd() * "/LES/" * case * "_profiles.jld2"
            # construct default loss function
            N = resolution[1]
            Δt = resolution[2]
            ℒ = closure_fixed_time_loss_function(filename; N=resolution[1], Δt = resolution[2], final_day = day)
            # choose default parameters
            optimal_𝑪 = copy(default_𝑪)
            # optimize using default optimize, iterate a few times
            println("-------------------")
            println("For case $case ")
            println("and resolution " * string(resolution[1]))
            println("and day $day")
            println("starting mcmc")
            Random.seed!(1234)
            # create name
            resolution_label = "_res_" * string(resolution[1])
            time_label = "_time_" * @sprintf("%.2f ", day)
            filename = pwd() * "/mcmc_data/" * case * resolution_label * time_label *  "_optima.jld2"
            # end of name creation
            mcmc_data = jldopen(filename, "r")
            initial_𝑪 = mcmc_data["parameter"]
            ℒ⁰ = mcmc_data["loss"]
            close(mcmc_data)
            # scale the loss function by ℒ
            nll(𝑪) = ℒ(𝑪) / ℒ⁰
            filename = pwd() * "/mcmc_data/" * case * resolution_label * time_label *  "_mcmc.jld2"
            # parameters for mcmc
            nt = 10000
            frequency = 100
            # define proposal matrix, 5% of default value
            proposal =  CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
            # now markov chain
            CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename, verbose = true)
            println("done")
        end
    end
end
