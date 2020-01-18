# Non nonlocal term comparison plot
include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
using Plots, Printf, Statistics, JLD2

# use PyPlot backend
pyplot()
# local vs nonlocal kpp figures

save_figures = true

# choose case
case = cases[1]

# get LES
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

# define loss functions and forward maps
subsample = 1:1:length(les.t)
N = 16
Δt = 10 * 60 #seconds
zᵖ = zeros(N)
# define the forward map
𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample, grid = zᵖ)
# define the loss function
ℒ = CoreFunctionality.closure_T_nll(𝒢, les)
# define time dependent loss function
ℒᵗ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )

# get MCMC data

filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
acceptance_rate = sum(e1 .== e2) / length(e1)
indmin = argmin(e1)
close(mcmc_data)
# get MCMC nonolocal data
resolution_label = "_res_" * string(N)
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_no_nonlocal_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
nn_chain = mcmc_data["𝑪"]
nn_e1 = mcmc_data["ε"]
nn_e2 = mcmc_data["proposal_ε"]
acceptance_rate = sum(e1 .== e2) / length(e1)
nn_indmin = argmin(nn_e1)
close(mcmc_data)

seconds_in_a_day = 86400
# parameters to loop over
p = []
labels = ["Default KPP", "Optimized KPP", "Local KPP"]
default_𝑪 = [0.1, 6.33, 1.36, 3.19]
optimal_𝑪 = chain[:, indmin]
nn_optimal_𝑪 = nn_chain[:, nn_indmin]
parameter_list =[default_𝑪, optimal_𝑪, nn_optimal_𝑪]
plot()
for j in 1:3
    𝑪 = parameter_list[j]
    loss = ℒ(𝑪)
    Tᵖ = 𝒢(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les.T[:,end], les.z, label = "LES", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p,p1)
end

p1 = plot(p[2:3]...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_5.pdf")
end

###
# loss as a function in time
p2 = plot()
loss_p = []
for j in 1:3
    𝑪 = parameter_list[j]
    t = les.t ./ seconds_in_a_day
    inds = 30:length(les.t)
    loss = ℒᵗ(𝑪)
    p2 = plot!(t[inds], sqrt.(loss[inds]), label = labels[j], legend = :topleft, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    push!(loss_p, p2)
end

plot(loss_p[2])

if save_figures == true
    savefig(p2, pwd() * "/figures/figure_5_alternate.pdf")
end


###
# create .gif
plot()
anim = @animate for i in 1:20:length(les.t)
    p = []
    for j in 1:3
        𝑪 = parameter_list[j]
        loss = ℒᵗ(𝑪)
        Tᵖ = 𝒢(𝑪)
        loss_string = @sprintf("%.1e", sqrt(loss[i]))
        p1 = plot(les.T[:,i], les.z, label = "LES", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
        p1 = scatter!(Tᵖ[:,i], zᵖ, label = labels[j], title = "Error = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
        push!(p,p1)
    end
    display(plot(p...))
end

if save_figures == true
    gif(anim, pwd() * "/figures/figure_5_dynamic.gif", fps = 15)
end
