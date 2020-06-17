include(pwd() * "/src/LocalOceanUQSupplementaryMaterials.jl")
include(pwd() * "/scripts/utils.jl")

generate_example = false
generate_plot = true

use_covariance_estimate = true
case_range = 1:2
resolution_range = 1:1
const ensemble_size = 10^5

resolution = resolutions[1]
case = cases[1]

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

### Run
if generate_example
    # to reset rightbounds
    include(pwd() * "/scripts/utils.jl")
    # scale the loss function by ℒ
    Random.seed!(12345)
    const factor = 10
    nll(𝑪) = ℒ(𝑪) / ℒ⁰ / factor
    filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
    # parameters for mcmc
    nt = ensemble_size
    frequency = 100
    # define proposal matrix, 5% of default value
    σ[1] = 10^(-6)
    σ[2] = 10^(-6)
    σ[3] = 0.1
    σ[4] = 0.3
    proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
    use_covariance_estimate = false
    if use_covariance_estimate
        proposal = CoreFunctionality.closure_proposal(Σ, left_bounds = left_bounds, right_bounds = right_bounds)
    end
    # now markov chain
    guess_𝑪 = copy(initial_𝑪)
    guess_𝑪[3] += 2.0
    guess_𝑪[4] += 6.0
    CoreFunctionality.markov_chain(nll, guess_𝑪, proposal, nt,  freq = frequency, filename = filename, verbose = true)
    println("done")
end
###
pyplot(size = (400,400))
if generate_plot
    const factor = 10
    filename = filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
    mcmc_data = jldopen(filename, "r")
    chain = mcmc_data["𝑪"]
    e1 = mcmc_data["ε"]
    e2 = mcmc_data["proposal_ε"]
    acceptance_rate = sum(e1 .== e2) / length(e1)
    println("the acceptance rate was")
    println(acceptance_rate)
    indmin = argmin(e1)
    close(mcmc_data)
    index1 = 3
    index2 = 4
    bools = e1 .< minimum(e1) * 2
    tmp_ind = argmax(bools)
    if factor > 1
        tmp_ind = 80
        # tmp_ind = 830
    else
        tmp_ind = 100
        # tmp_ind = 1037
    end
    bins = 200
    Cᴿ = 0.3
    chain[4, :] *= Cᴿ
    right_bounds[4] = 4.0

    p1 = histogram2d(chain[index1, tmp_ind:end], chain[index2, tmp_ind:end], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, legend = false, color = cgrad(:blues, rev = false), textsize = 10)
    # another way to accomplish similar things is with
    # marker_z = (+), color = :bluesreds
    # see http://docs.juliaplots.org/latest/generated/plotly/#plotly-ref35-1
    for i in 1:tmp_ind
        ω = i / tmp_ind / 8 * 4.5
        p1 = scatter!(chain[index1, i:i], chain[index2, i:i], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 1.0, RGB(0.1,ω,1-ω), stroke(1, 1.0, :black, :dot)), label = false)
    end
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = "starting value", markersize= 15)
    scatter!(initial_𝑪[index1, 1:1], initial_𝑪[index2, 1:1] .* Cᴿ, shape = :star, color = :green, label = "optimal value", legend = :topright, markersize= 15, legendfont = font("Times new roman", 13))
    display(p1)
    savefig(p1, pwd() * "/figures/simpler_mcmc_"* string(factor) * "_v2.pdf")
end

if generate_plot
    index = 3
    Δx = right_bounds[index] - left_bounds[index]
    Δy = length(chain[index, tmp_ind:end-1]) / 10
    ratio = 1/3 * Δx / Δy
    p2 = histogram(chain[index, tmp_ind:end-1],  bins = bins, legend = false, normalize = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :none, aspect_ratio = ratio, xlims = (left_bounds[index], right_bounds[index]), ylims = (0, Δy), ticks = false, edges = false, linewidth = 0.1)
    display(p2)
    savefig(p2, pwd() * "/figures/simpler_mcmc_"*string(factor)*"_marginal_"*string(index)*".pdf")
end

if generate_plot
    index = 4
    Δx = right_bounds[index] - left_bounds[index]
    Δy = 2.0
    ratio = 1/3 * Δx / Δy
    p2 = histogram(chain[index, tmp_ind:end-1],  bins = bins, legend = false, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :none, aspect_ratio = ratio, xlims = (left_bounds[index], right_bounds[index]), ylims = (0, Δy), ticks = false, linewidth = 0.1)
    savefig(p2, pwd() * "/figures/simpler_mcmc_" * string(factor) * "_marginal_"*string(index)*".pdf")
    display(p2)
end

###
pp1 = histogram(chain[index, tmp_ind:end-1],  bins = 100, legend = false, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index], right_bounds[index]), ticks = true, linewidth = 0.1)

plot(pp1, pp2)
