include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")
using Plots, Printf, Statistics, JLD2, MCMCDiagnostics
# use PyPlot backend
pyplot()
#optimized vs nonoptimized kpp figures

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
resolution_label = "_res_" * string(N)
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
acceptance_rate = sum(e1 .== e2) / length(e1)
println("the acceptance rate was")
println(acceptance_rate)
indmin = argmin(e1)
close(mcmc_data)
ess = randn(4)
for i in 1:4
    # one million is a bit much
    x1 = chain[i,1:end]
    variance_x1 = var(x1)
    ess[i] = effective_sample_size(x1, variance_x1)
end
println("the effective sample size was")
println(ess)

seconds_in_a_day = 86400
# parameters to loop over
labels = ["Default", "Mode", "Mean", "Median"]
default_𝑪 = [0.1, 6.33, 1.36, 3.19]
#default_𝑪 = [0.0760809666611145; 4.342473912404762; 2.1630355831002954; 5.57111619953263] # from megachain
optimal_𝑪 = chain[:, indmin]
mean_𝑪 = mean(chain,dims=2)
median_𝑪 = median(chain,dims=2)
#median_𝑪 = [0.0760809666611145; 4.342473912404762; 2.1630355831002954; 5.57111619953263] # across all

parameter_list =[default_𝑪, optimal_𝑪, mean_𝑪, median_𝑪]
plot()
p = []
for j in 1:4
    𝑪 = parameter_list[j]
    loss = ℒ(𝑪)
    Tᵖ = 𝒢(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les.T[:,end], les.z, label = "LES", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius)
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error = " * loss_string * " " * celsius, markersize = 3, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p,p1)
end

p1 = plot(p[1:2]...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_2.png")
end
p1 = plot(p...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_2b.png")
end

###
# loss as a function in time
p2 = plot()
loss_p = []
for j in 1:4
    𝑪 = parameter_list[j]
    t = les.t ./ seconds_in_a_day
    inds = 30:length(les.t)
    loss = ℒᵗ(𝑪)
    p2 = plot!(t[inds], sqrt.(loss[inds]), label = labels[j], legend = :topleft, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    push!(loss_p, p2)
end

plot(loss_p[2])

if save_figures == true
    savefig(p2, pwd() * "/figures/figure_2_alternate.png")
end


###
# create .gif
plot()
anim = @animate for i in 1:20:length(les.t)
    p = []
    for j in 1:4
        𝑪 = parameter_list[j]
        loss = ℒᵗ(𝑪)
        Tᵖ = 𝒢(𝑪)
        loss_string = @sprintf("%.1e", sqrt(loss[i]))
        p1 = plot(les.T[:,i], les.z, label = "LES", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius)
        p1 = scatter!(Tᵖ[:,i], zᵖ, label = labels[j], title = "Error = " * loss_string * " " * celsius, markersize = 3, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
        push!(p,p1)
    end
    display(plot(p...))
end
if save_figures == true
    gif(anim, pwd() * "/figures/figure_2_dynamic.gif", fps = 15)
end
