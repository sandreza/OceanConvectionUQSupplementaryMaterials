include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")

# use optimal parameters as starting point for next itration
parameter_continuation = true
random_starting_point_1 = false
random_starting_point_2 = false
σ = default_𝑪 * 0.5 # more liberal searching for optimal parameters
case_range = 1:1
for resolution in resolutions[1:1]
    # reset default parameters if parameter continuation is used
    default_𝑪 = [0.1, 6.33, 8.36, 3.19]
    # default_𝑪 = [0.11803164331592443, 3.7246545857676954, 0.35191154207167974, 6.225750233165317] #best for 1
    # default_𝑪 = [0.04874744540063653, 3.760819427517219, 0.1814772890705244, 11.98844037974979]   #best for 2
    σ = default_𝑪 * 0.5
    if random_starting_point_1
        default_𝑪[1] = rand()
        default_𝑪[2] = rand() * 8
        default_𝑪[3] = rand() * 6
        default_𝑪[4] = rand() * 16
    end
    for case in cases[case_range]
        if random_starting_point_2
            default_𝑪[1] = rand()
            default_𝑪[2] = rand() * 8
            default_𝑪[3] = rand() * 6
            default_𝑪[4] = rand() * 16
        end
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        # construct default loss function
        N = resolution[1]
        Δt = resolution[2]
        nll = closure_default_loss_function(filename, N = N, Δt = Δt)
        # optimize using default optimize, iterate a few times
        println("-------------------")
        println("For case $case ")
        println("And KPP gridpoints N = $N")
        println("starting optimization with default parameters")
        println(default_𝑪)
        proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
        # for reproducibility set random seed
        Random.seed!(1234)
        default_ℒ = nll(default_𝑪)
        # random walk optimization
        println("random walk")
        # optimal_𝑪 = CoreFunctionality.optimize(default_𝑪, nll; nt = 1000, restart = 2, proposal = proposal, rescale = true, freq = 100, scale = default_ℒ)
        optimal_𝑪, Σ = CoreFunctionality.optimize_and_estimate_proposal(default_𝑪, nll, left_bounds, right_bounds, nt = 2000, restart = 3, proposal = [], filename = [], rescale = true, freq = 100, verbose = true)

        default_ℒ = nll(default_𝑪)
        optimal_ℒ = nll(optimal_𝑪)
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
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_optima.jld2"
        parameter = optimal_𝑪
        loss = optimal_ℒ
        covariance = Σ
        @save filename parameter loss covariance
        # parameter continuation
        if parameter_continuation == true
            # constant extrapolation
            @. default_𝑪 = optimal_𝑪
            # linear extrapolation, useful when varying stratificaton first
            # @. default_𝑪 = 2 * optimal_𝑪 - default_𝑪
        end
    end
end
