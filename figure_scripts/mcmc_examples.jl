include(pwd() * "/src/LocalOceanUQSupplementaryMaterials.jl")
include(pwd() * "/scripts/utils.jl")

generate_example = false
generate_plot = true

use_covariance_estimate = true
case_range = 1:2
resolution_range = 1:1


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
    const ensemble_size = 10^5
    inverse_factor = true
    include(pwd() * "/scripts/utils.jl")
    # scale the loss function by ℒ
    Random.seed!(12345)
    const factor = 10
    if inverse_factor
        nll(𝑪) = ℒ(𝑪) / ℒ⁰ * factor
    else
        nll(𝑪) = ℒ(𝑪) / ℒ⁰ / factor
    end
    if inverse_factor
    filename = pwd() * "/mcmc_data/" * "toy_example_i" *  string(factor) * "_mcmc.jld2"
    else
        filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
    end
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
# Version 1
pyplot(size = (500, 500))
if generate_plot
    const factor = 10
    inverse_factor = true
    if inverse_factor
        filename = filename = pwd() * "/mcmc_data/" * "toy_example_i" *  string(factor) * "_mcmc.jld2"
    else
        filename = filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
    end
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
        tmp_ind = 830
    else
        tmp_ind = 100
        tmp_ind = 1037
    end
    if inverse_factor
        tmp_ind = 1000
    end
    bins = 200
    Cᴿ = 0.3
    chain[4, :] *= Cᴿ
    right_bounds[4] = 4.0

    p1 = histogram2d(chain[index1, tmp_ind:end], chain[index2, tmp_ind:end], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, legend = false, color = cgrad(:blues, rev = false), xtickfont=font(18), ytickfont=font(18), xguidefontsize=18, yguidefontsize = 18, legendfontsize = 18)
    # Starting value
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = "starting value", markersize= 15, opacity = 0.5)
    # another way to accomplish similar things is with
    # marker_z = (+), color = :bluesreds
    # see http://docs.juliaplots.org/latest/generated/plotly/#plotly-ref35-1
    for i in 1:tmp_ind
        ω = i / tmp_ind / 8 * 4.5
        p1 = scatter!(chain[index1, i:i], chain[index2, i:i], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 0.5, RGB(0.1,ω,1-ω), stroke(0.1, 0.1, :black, :dot)), label = false)
    end
    i = tmp_ind +1
    ω = i / tmp_ind / 8 * 4.5
    p1 = scatter!(chain[index1, i:i], chain[index2, i:i], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 0.5, RGB(0.1,ω,1-ω), stroke(1, 1.0, :black, :dot)), label = "chain")
    scatter!(initial_𝑪[index1, 1:1], initial_𝑪[index2, 1:1] .* Cᴿ, shape = :star, color = :green, label = "optimal value", legend = :topright, markersize= 15, legendfont = font("Times new roman", 13), opacity = 0.5)
    display(p1)
    if inverse_factor
        savefig(p1, pwd() * "/figures/simpler_mcmc_i"* string(factor) * "_v1.pdf")
    else
        savefig(p1, pwd() * "/figures/simpler_mcmc_"* string(factor) * "_v1.pdf")
    end
end

### Another version of the plots
# Version 2
pyplot(size = (500,500))
if generate_plot
    const factor = 10
    inverse_factor = true
    if inverse_factor
        filename = filename = pwd() * "/mcmc_data/" * "toy_example_i" *  string(factor) * "_mcmc.jld2"
    else
        filename = filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
    end
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
    if inverse_factor
        tmp_ind = 1000
    end
    bins = 300
    Cᴿ = 0.3
    chain[4, :] *= Cᴿ
    right_bounds[4] = 4.0

    p1 = histogram2d(chain[index1, tmp_ind:end], chain[index2, tmp_ind:end], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, legend = false, color = cgrad(:blues, rev = false), textsize = 10, xtickfont=font(18), ytickfont=font(18), xguidefontsize=18, yguidefontsize = 18, legendfontsize = 18)
    # Starting Value
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = "starting value", markersize= 15, opacity = 0.5)
    # Burn-in
    p1 = scatter!(chain[index1, 1:tmp_ind], chain[index2, 1:tmp_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 0.25, RGB(1.0, 0.30, 0.0), stroke(1, 1.0, :black, :dot)), label = "burn-in")
    # To plot over burn in
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = false, markersize= 15, opacity = 0.5)
    # Best-known optimal value
    scatter!(initial_𝑪[index1, 1:1], initial_𝑪[index2, 1:1] .* Cᴿ, shape = :star, color = :green, label = "optimal value", legend = :topright, markersize= 15, legendfont = font("Times new roman", 13), opacity = 0.5)
    display(p1)
    if inverse_factor
        savefig(p1, pwd() * "/figures/simpler_mcmc_i"* string(factor) * "_v2.pdf")
    else
        savefig(p1, pwd() * "/figures/simpler_mcmc_"* string(factor) * "_v2.pdf")
    end
end

###
# Version 3
pyplot(size = (500,500))
if generate_plot
    const factor = 10
    inverse_factor = false
    if inverse_factor
        filename = filename = pwd() * "/mcmc_data/" * "toy_example_i" *  string(factor) * "_mcmc.jld2"
    else
        filename = filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
    end
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
        tmp_ind = 80 # burn-in
        tmp_ind = 830
    else
        tmp_ind = 100 # burn-in
        tmp_ind = 1037
    end
    if inverse_factor
        tmp_ind = 1100
    end
    bins = 300
    Cᴿ = 0.3
    chain[4, :] *= Cᴿ
    right_bounds[4] = 4.0

    p1 = histogram2d(chain[index1, tmp_ind:end], chain[index2, tmp_ind:end], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, legend = false, color = cgrad(:blues, rev = false), textsize = 10, xtickfont=font(18), ytickfont=font(18), xguidefontsize=18, yguidefontsize = 18, legendfontsize = 18)
    # Starting Value
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = "starting value", markersize= 15, opacity = 0.5)
    # Burn-in
    p1 = scatter!(chain[index1, 1:1:tmp_ind], chain[index2, 1:1:tmp_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 0.25, RGB(1.0, 0.30, 0.0), stroke(1, 1.0, :black, :dot)), label = "chain")
    # Best-known optimal value
    scatter!(initial_𝑪[index1, 1:1], initial_𝑪[index2, 1:1] .* Cᴿ, shape = :star, color = :green, label = "optimal value", legend = :topright, markersize= 15, legendfont = font("Times new roman", 13), opacity = 0.5)
    display(p1)
    if inverse_factor
        savefig(p1, pwd() * "/figures/simpler_mcmc_i"* string(factor) * "_v3.pdf")
    else
        savefig(p1, pwd() * "/figures/simpler_mcmc_"* string(factor) * "_v3.pdf")
    end
end


###
const factor = 10
inverse_factor = true
if inverse_factor
    filename = filename = pwd() * "/mcmc_data/" * "toy_example_i" *  string(factor) * "_mcmc.jld2"
else
    filename = filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
end
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
acceptance_rate = sum(e1 .== e2) / length(e1)
println("the acceptance rate was")
println(acceptance_rate)
indmin = argmin(e1)
close(mcmc_data)
Cᴿ = 0.3
chain[4, :] *= Cᴿ
right_bounds[4] = 4.0

# Marginals
if generate_plot
    bins = 60
    index = 3
    Δx = right_bounds[index] - left_bounds[index]
    Δy = 7
    ratio = 1/3 * Δx / Δy
    p2 = histogram(chain[index, tmp_ind:end-1],  bins = bins, legend = false, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, xlims = (left_bounds[index], right_bounds[index]), ylims = (0, Δy), ticks = false, edges = false, linewidth = 0.00)
    display(p2)
    if inverse_factor
        savefig(p2, pwd() * "/figures/simpler_mcmc_"*string(factor)*"_marginal_i"*string(index)*".pdf")
    else
        savefig(p2, pwd() * "/figures/simpler_mcmc_"*string(factor)*"_marginal_"*string(index)*".pdf")
    end
end

if generate_plot
    index = 4
    Δx = right_bounds[index] - left_bounds[index]
    Δy = 5
    ratio = 1/3 * Δx / Δy
    p2 = histogram(chain[index, tmp_ind:end-1],  bins = bins, legend = false, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = ratio, xlims = (left_bounds[index], right_bounds[index]), ylims = (0, Δy), ticks = false, linewidth = 0.0)
    if inverse_factor
        savefig(p2, pwd() * "/figures/simpler_mcmc_"*string(factor)*"_marginal_i"*string(index)*".pdf")
    else
        savefig(p2, pwd() * "/figures/simpler_mcmc_"*string(factor)*"_marginal_"*string(index)*".pdf")
    end
    display(p2)
end


###
# gif for fun, takes a little while to run
const factor = 10
inverse_factor = false
if inverse_factor
    filename = filename = pwd() * "/mcmc_data/" * "toy_example_i" *  string(factor) * "_mcmc.jld2"
else
    filename = filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
end
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
proposal_chain = mcmc_data["proposal_𝑪"]
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
    tmp_ind = 830
else
    tmp_ind = 100
    tmp_ind = 1037
end
if inverse_factor
    tmp_ind = 1000
end
bins = 200
Cᴿ = 0.3
chain[4, :] *= Cᴿ
proposal_chain[4, :] *= Cᴿ
right_bounds[4] = 4.0
using LinearAlgebra
tail_ind = 10
anim = @animate for i in 1:1:2*10^3
    p1 = histogram2d(chain[index1, tmp_ind:end], chain[index2, tmp_ind:end], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, legend = false, color = cgrad(:blues, rev = false), xtickfont=font(18), ytickfont=font(18), xguidefontsize=18, yguidefontsize = 18, legendfontsize = 18)
    # Starting value
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = "starting value", markersize= 15, opacity = 0.5)
    # another way to accomplish similar things is with
    # marker_z = (+), color = :bluesreds
    # see http://docs.juliaplots.org/latest/generated/plotly/#plotly-ref35-1
    ω = i / tmp_ind / 8 * 4.5
    p1 = scatter!(chain[index1, i:i+tail_ind], chain[index2, i:i+tail_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 0.5, RGB(0.0, 0.0, 1.0), stroke(0.1, 0.1, :black, :dot)), label = false)

    scatter!(chain[index1, i+1+tail_ind:i+1+tail_ind], chain[index2, i+1+tail_ind:i+1+tail_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 1.0, RGB(0.0, 0.0, 1.0), stroke(0.1, 0.1, :black, :dot)), label = "chain")
    # proposal parameter
    accepted = norm(proposal_chain[:, i+2+tail_ind] - chain[:, i+2+tail_ind]) < eps(1.0)
    if accepted
        color = RGB(0.0, 1.0, 0.0)
    else
        color = RGB(1.0, 0.0, 0.0)
    end
    scatter!(proposal_chain[index1, i+2+tail_ind:i+2+tail_ind], proposal_chain[index2, i+2+tail_ind:i+2+tail_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 1.0, color, stroke(1, 1.0, :black, :dot)), label = "proposal")
    # end location
    scatter!(initial_𝑪[index1, 1:1], initial_𝑪[index2, 1:1] .* Cᴿ, shape = :star, color = :green, label = "optimal value", legend = :topright, markersize= 15, legendfont = font("Times new roman", 13), opacity = 0.5)
end

if inverse_factor
    gif(anim, pwd() * "/" * string(factor) * "_i_scratch.gif", fps = 15)
else
    gif(anim, pwd() * "/" * string(factor) * "_scratch.gif", fps = 15)
end
