include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")

# boolean label for optimizing or performing mcmc
optimize_compromise = false
mcmc_compromise = false

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
    optimal_𝑪 = ( 𝑪1 .+ 𝑪2 ) ./ 2
    println("starting optimization")
    for j in 1:5
        println("loop " * string(j))
        res = optimize(ℒ_compromise, optimal_𝑪)
        tmp_𝑪  = Optim.minimizer(res)
        @. optimal_𝑪 = tmp_𝑪
    end
    default_𝑪  = ( 𝑪1 .+ 𝑪2 ) ./ 2
    res = optimize(ℒ_compromise, optimal_𝑪)
    tmp_𝑪  = Optim.minimizer(res)
    optimal_𝑪 = Optim.minimizer(res)
    optimal_ℒ = Optim.minimum(res)
    default_ℒ = ℒ_compromise(default_𝑪)
    #print info
    println("The default parameters are $default_𝑪")
    println("The lossfunction value is $(default_ℒ)")
    println("The optimal parameters are $optimal_𝑪")
    println("The lossfunction value is $optimal_ℒ")
    println("The improvement is a factor of $(default_ℒ/optimal_ℒ)")
    println("-------------------")

    # save optimal values and loss function value
    resolution_label = "_res_" * string(resolution[1])
    filename = pwd() * "/mcmc_data/" * case_name * resolution_label * "_optima.jld2"
    parameter = optimal_𝑪
    loss = optimal_ℒ
    @save filename parameter loss
end

if mcmc_compromise == true
    filename = pwd() * "/mcmc_data/" * case_name * resolution_label * "_optima.jld2"
    mcmc_data = jldopen(filename, "r")
    initial_𝑪 = mcmc_data["parameter"]
    ℒ⁰ = mcmc_data["loss"]
    close(mcmc_data)
    # scale the loss function by ℒ
    nll(𝑪) = ℒ_compromise(𝑪) / ℒ⁰
    filename = pwd() * "/mcmc_data/" * case_name * resolution_label * "_mcmc.jld2"
    # parameters for mcmc
    nt = 20000
    frequency = 100
    # define proposal matrix, 5% of default value
    proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
    # now markov chain
    CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename)
    println("done")
end
