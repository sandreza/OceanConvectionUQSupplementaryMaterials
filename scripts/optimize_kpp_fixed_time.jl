include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")

# use optimal parameters as starting point for next itration
parameter_continuation = true

case_range = 1:1
resolution_range = 1:1
days = [1.0, 2.0 , 3.0 , 4.0 , 5.0, 6.0, 7.0, 8.0]

for resolution in resolutions[resolution_range]
    # reset default parameters if parameter continuation is used
    default_𝑪 = [0.1, 6.33, 1.36, 3.19]
    for case in cases[case_range]
        for day in days
            filename = pwd() * "/LES/" * case * "_profiles.jld2"
            # construct default loss function
            N = resolution[1]
            Δt = resolution[2]
            nll = closure_fixed_time_loss_function(filename; N=resolution[1], Δt = resolution[2], final_day = day)
            # optimize using default optimize, iterate a few times
            println("-------------------")
            println("For case $case ")
            println("And KPP gridpoints N = $N")
            println("And day $day")
            println("starting optimization with default parameters")
            println(default_𝑪)
            proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
            # for reproducibility set random seed
            Random.seed!(1234)
            default_ℒ = nll(default_𝑪)
            # random walk optimization
            println("random walk")
            optimal_𝑪 = CoreFunctionality.optimize(default_𝑪, nll; nt = 1000, restart = 2, proposal = proposal, rescale = true, freq = 100, scale = default_ℒ)

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
            time_label = "_time_" * @sprintf("%.2f ", day)
            filename = pwd() * "/mcmc_data/" * case * resolution_label * time_label *  "_optima.jld2"
            # done with filename create
            parameter = optimal_𝑪
            loss = optimal_ℒ
            @save filename parameter loss
            # parameter continuation
            if parameter_continuation == true
                # constant extrapolation
                @. default_𝑪 = optimal_𝑪
                # linear extrapolation, useful when varying stratificaton first
                # @. default_𝑪 = 2 * optimal_𝑪 - default_𝑪
            end
        end
    end
end
