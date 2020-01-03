include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")

# use optimal parameters as starting point for next itration
parameter_continuation = false
σ = default_𝑪 * 0.2 # more liberal searching for optimal parameters
case_range = 20:1:34
case_range = [1,6,7,8,9,10]
case_range = 2:2
#case_range = 7:10
for resolution in resolutions[1:1]
    for case in cases[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        # construct default loss function
        les = CoreFunctionality.OceananigansData(filename)

        N = resolution[1]
        Δt = resolution[2]
        nll = closure_flexible_loss_function(filename, N = N, Δt = Δt, power = 1.0)
        # optimize using default optimize, iterate a few times
        println("-------------------")
        println("For case $case ")
        println("And KPP gridpoints N = $N")
        println("starting optimization with default parameters")
        resolution_label = "_res_" * string(N)

        NN = sqrt(les.α * les. g * les.bottom_T)
        default_𝑪 = [1e-4, 3.5 * 1.0, 10.0, 0.0, 0.375, NN]
        default_𝑪 = [0.1, 3.5 * 1.0, 1.0, 0.0, 0.375, NN]
        default_𝑪 = [0.005741998337334633, 3.629207116893695, 1.1392751590144323, 0.0, 0.40974149273298843, NN]
        println(default_𝑪)
        σ = default_𝑪 * 0.5
        σ[6] = eps(1.0)
        #=
        left_bounds = [0.0, 3.0, 5.0, 0.0, 0.0, NN]
        right_bounds = [0.01, 5.0, 10.0, eps(1.0), 1.0, NN + eps(1.0)]
        =#
        left_bounds = [0.0, 0.0, 0.0, 0.0, 0.0, NN]
        right_bounds = [0.01, 8.0, 10.0, eps(1.0), 1.0, NN + eps(1.0)]

        proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
        # for reproducibility set random seed
        Random.seed!(1234)
        default_ℒ = nll(default_𝑪)
        # random walk optimization
        println("random walk")
        # optimal_𝑪 = CoreFunctionality.optimize(default_𝑪, nll; nt = 1000, restart = 0, proposal = proposal, rescale = true, freq = 100, scale = default_ℒ)
        optimal_𝑪, Σ = CoreFunctionality.optimize_and_estimate_proposal(default_𝑪, nll, left_bounds, right_bounds, nt = 1000, restart = 1, proposal = [], filename = [], rescale = true, freq = 100, verbose = true)

        default_ℒ = nll(default_𝑪)
        optimal_ℒ = nll(optimal_𝑪)
        println("--------------------")
        println("The default parameters are $default_𝑪")
        println("The lossfunction value is $(default_ℒ)")
        println("The optimal parameters are $optimal_𝑪")
        println("The lossfunction value is $optimal_ℒ")
        println("The improvement is a factor of $(default_ℒ/optimal_ℒ)")
        println("-------------------")
        # save optimal values and loss function value
        resolution_label = "_res_" * string(resolution[1])
        extra_label = "_flexible_new"
        filename = pwd() * "/mcmc_data/" * case * resolution_label * extra_label *  "_optima.jld2"
        parameter = optimal_𝑪
        loss = optimal_ℒ
        covariance = Σ
        @save filename parameter loss covariance
    end
end
