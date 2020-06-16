# for John
include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

generate_example = false
generate_plot = true

generate_simple_example = false
generate_simple_plot = true
use_covariance_estimate = true

case = cases[6]
if generate_example
    resolution = resolutions[1]
    use_covariance_estimate = true
    # construct filename
    filename = pwd() * "/LES/" * case * "_profiles.jld2"
    # construct default loss function
    N = resolution[1]
    Δt = resolution[2]
    ℒ = closure_default_loss_function(filename, N = N, Δt = Δt)
    # optimize using default optimize, iterate a few times
    println("-------------------")
    println("For case $case ")
    println("and resolution " * string(resolution[1]))
    println("starting mcmc")
    Random.seed!(1234)
    resolution_label = "_res_" * string(resolution[1])
    filename = pwd() * "/mcmc_data/" * case * resolution_label * "_optima.jld2"
    mcmc_data = jldopen(filename, "r")
    initial_𝑪 = mcmc_data["parameter"] .+ [0.0 , 0.0, 1.0, 4.0]
    # [0.2306826206143779, 3.8047948889896963, 3.0, 12.0]
    ℒ⁰ = mcmc_data["loss"]
    if use_covariance_estimate
        Σ = mcmc_data["covariance"]
    end
    close(mcmc_data)
    # scale the loss function by ℒ
    nll(𝑪) = ℒ(𝑪) / ℒ⁰
    filename = pwd() * "/mcmc_data/" * case * resolution_label * "_example_mcmc.jld2"
    # parameters for mcmc
    nt = 1000
    frequency = 100
    # define proposal matrix, 5% of default value
    proposal = CoreFunctionality.closure_proposal(σ, left_bounds = left_bounds, right_bounds = right_bounds)
    if use_covariance_estimate
        Σ *= 0.01
        proposal = CoreFunctionality.closure_proposal(Σ, left_bounds = left_bounds, right_bounds = right_bounds)
        println("using covariance estimate")
    end
    # now markov chain
    CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename, verbose = true)
    println("done")
end

###
resolution = resolutions[1]
resolution_label = "_res_" * string(resolution[1])
if generate_plot
    filename = pwd() * "/mcmc_data/" * case * resolution_label * "_example_mcmc.jld2"
    mcmc_data = jldopen(filename, "r")
    chain = mcmc_data["𝑪"]
    e1 = mcmc_data["ε"]
    e2 = mcmc_data["proposal_ε"]
    acceptance_rate = sum(e1 .== e2) / length(e1)
    println("the acceptance rate was")
    println(acceptance_rate)
    indmin = argmin(e1)
    close(mcmc_data)
    gr()
    pair = [3, 4]

    index1 = pair[1]
    index2 = pair[2]
    bins = 30
    p1 = histogram2d(chain[index1, :], chain[index2, :], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :red, markersize = 6, label = "start")
    # scatter!(chain[index1, end:end], chain[index2, end:end], shape = :star, color = :yellow, markersize = 8, label = "optimal")
    display(p1)
    savefig(p1, pwd() * "/figures/kpp_mcmc.pdf")
end


###
# simpler example
Random.seed!(1234)
const μ_exact = [1.0; 2.0]
Σ = [1 1/2; 1/2 1] .* 1.0
iΣ = inv(Σ) .* 0.5 # for simplicity put the 1/2 factor here
pdf(x) = exp( - (x-μ_exact)' * iΣ * (x-μ_exact) )

initial_𝑪 = [3.0; 6.0]
nll(𝑪) = - log(pdf(𝑪))

if generate_simple_example
    initial_𝑪 = [3.0; 6.0] * 4.0
    nll(𝑪) = - log(pdf(𝑪))
    # parameters for mcmc
    nt = 100000
    frequency = 100
    # define proposal matrix, 5% of default value
    Σ_p = [0.02 0.0; 0.0 0.02]
    proposal = CoreFunctionality.closure_proposal(Σ_p)
    # now markov chain
    filename = pwd() * "/mcmc_data/"* "simple_example_mcmc.jld2"
    CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename, verbose = true)
end
###
if generate_simple_plot
    filename = pwd() * "/mcmc_data/"* "simple_example_mcmc.jld2"
    mcmc_data = jldopen(filename, "r")
    chain = mcmc_data["𝑪"]
    e1 = mcmc_data["ε"]
    e2 = mcmc_data["proposal_ε"]
    acceptance_rate = sum(e1 .== e2) / length(e1)
    println("the acceptance rate was")
    println(acceptance_rate)
    indmin = argmin(e1)
    close(mcmc_data)

    pyplot()
    pair = [1, 2]
    index1 = pair[1]
    index2 = pair[2]
    bins = 250
    p1 = histogram2d(chain[index1, :], chain[index2, :], xlabel = "parameter 1", ylabel = "parameter 2", bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :star, color = :blue, label = "starting value", markersize= 15)
    scatter!(μ_exact[1:1], μ_exact[2:2], shape = :star, color = :green, label = "optimal value", legend = :topleft, markersize= 15, legendfont = font("Times new roman", 13))
    display(p1)
    savefig(p1, pwd() * "/figures/simple_mcmc.pdf")
end

###
# simpler example 2
Random.seed!(1234)
const μ_exact = [1.0; 2.0]
const example_scale = 1
Σ = [1 1/2; 1/2 1] .* example_scale
iΣ = inv(Σ) .* 0.5 # for simplicity put the 1/2 factor here
pdf(x) = exp( - (x-μ_exact)' * iΣ * (x-μ_exact) )

initial_𝑪 = [3.0; 6.0]
nll(𝑪) = - log(pdf(𝑪))

if generate_simple_example
    initial_𝑪 = [3.0; 6.0] * 4.0
    nll(𝑪) = - log(pdf(𝑪))
    # parameters for mcmc
    nt = 10^5
    frequency = 100
    # define proposal matrix, 5% of default value
    Σ_p = [0.02 0.0; 0.0 0.02] .* example_scale
    proposal = CoreFunctionality.closure_proposal(Σ_p)
    # now markov chain
    filename = pwd() * "/mcmc_data/"* "simple_example_mcmc_"*string(example_scale)*".jld2"
    CoreFunctionality.markov_chain(nll, initial_𝑪, proposal, nt,  freq = frequency, filename = filename, verbose = true)
end
###
pyplot(size = (500,500))
if generate_simple_plot
    filename = pwd() * "/mcmc_data/"* "simple_example_mcmc_"*string(example_scale)*".jld2"
    mcmc_data = jldopen(filename, "r")
    chain = mcmc_data["𝑪"]
    e1 = mcmc_data["ε"]
    e2 = mcmc_data["proposal_ε"]
    acceptance_rate = sum(e1 .== e2) / length(e1)
    println("the acceptance rate was")
    println(acceptance_rate)
    indmin = argmin(e1)
    close(mcmc_data)
    pair = [1, 2]
    index1 = pair[1]
    index2 = pair[2]
    bins = 250
    bools = e1 .< minimum(e1) * 100
    tmp_ind = argmax(bools)
    if example_scale > 1
        tmp_ind = 326
    else
        tmp_ind = 725
    end
    p1 = histogram2d(chain[index1, tmp_ind:end], chain[index2, tmp_ind:end], xlabel = "parameter 1", ylabel = "parameter 2", bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (-12, 13), ylims = (-10,25), color = cgrad(:blues, rev = false))
    # another way to accomplish similar things is with
    # marker_z = (+), color = :bluesreds
    # see http://docs.juliaplots.org/latest/generated/plotly/#plotly-ref35-1
    for i in 1:tmp_ind
        ω = i / tmp_ind / 8 * 4.5
        p1 = scatter!(chain[index1, i:i], chain[index2, i:i], xlabel = "parameter 1", ylabel = "parameter 2", bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (-12, 13), ylims = (-10,25), marker = (:hexagon, 4, 1.0, RGB(0.1,ω,1-ω), stroke(1, 1.0, :black, :dot)), label = false)
    end
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = "starting value", markersize= 15)
    scatter!(μ_exact[1:1], μ_exact[2:2], shape = :star, color = :green, label = "optimal value", legend = :topleft, markersize= 15, legendfont = font("Times new roman", 13))
    display(p1)
    savefig(p1, pwd() * "/figures/simple_mcmc_"*string(example_scale)*".pdf")
end

if generate_simple_plot
    index = 1
    a = [25, 35]
    b = [(-12, 13), (-10,25)]
    Δx = a[index]
    Δy = length(chain[index, 1:end-1]) / 20
    ratio = 1/3 * Δx / Δy
    p2 = histogram(chain[index, 1:end-1],  bins = bins, legend = false, normalize = false, ylabel = "arbitrary", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, xlims = b[index], ylims = (0, Δy), ticks = false)
    savefig(p2, pwd() * "/figures/simple_mcmc_"*string(example_scale)*"_marginal_"*string(index)*".pdf")
end

if generate_simple_plot
    index = 2
    a = [25, 35]
    b = [(-12, 13), (-10,25)]
    Δx = a[index]
    Δy = length(chain[index, 1:end-1]) / 20
    ratio = 1/3 * Δx / Δy
    p2 = histogram(chain[index, 1:end-1],  bins = bins, legend = false, normalize = false, ylabel = "arbitrary", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, xlims = b[index], ylims = (0, Δy), ticks = false)
    savefig(p2, pwd() * "/figures/simple_mcmc_"*string(example_scale)*"_marginal_"*string(index)*".pdf")
end
