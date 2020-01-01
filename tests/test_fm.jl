include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
using Plots, Printf, JLD2

animation = false
optimal = true
mega_chain_median = false
# Get LES data
wd = pwd()
case = cases[34]
filename = wd * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

# Set parameters
N = 16
Δt = 10*60; # 10 minutes

# define the forward map
zᵖ = zeros(N)
start = argmin(abs.(les.t .- 86400 * 0.25)) #start at a quarter of a day
#start = 1
subsample = start:1:length(les.t)
𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les,
                                              subsample = subsample, grid = zᵖ)

tˢ = collect(subsample) #indices for simulation time
finish = floor(Int, length(tˢ) / 1)

# compute the forward map
𝑪 = [0.1, 4.33, 1.36, 3.19 ] #/ 3.19 * 4.32]
#𝑪 = [0.48944664115621117, 5.606733551864946, 0.06980180915616958, 11.131016010149851]
if optimal
    resolution_label = "_res_" * string(N)
    filename2 = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
    mcmc_data = jldopen(filename2, "r")
    chain = mcmc_data["𝑪"]
    e1 = mcmc_data["ε"]
    e2 = mcmc_data["proposal_ε"]
    close(mcmc_data)
    𝑪 = chain[:, argmin(e1)]
    println(𝑪)
end

kpp_T = 𝒢(𝑪)


minT = minimum(les.T⁰)
maxT = maximum(les.T⁰)


if animation
    for j in 1:50:finish
        days = @sprintf("%.1f", les.t[j]/86400)
        jˢ = tˢ[j]
        Tˢ = les.T[:,jˢ]
        Tᵖ = kpp_T[:,j]
        z  = les.z

        p1 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT))
        scatter!( Tᵖ , zᵖ, label = "KPP", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
        display(p1)
    end
end

j = finish
jˢ = tˢ[j]
days = @sprintf("%.1f", les.t[jˢ]/86400)
Tˢ = les.T[:,jˢ]
Tᵖ = kpp_T[:,j]
z  = les.z

p1 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT))
scatter!( Tᵖ , zᵖ, label = "KPP", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
display(p1)

###
# Test other forward map

# Get LES data
#=
wd = pwd()
#3:8:34, 1, 9, 17, 25
case = cases2[1]
filename = wd * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)
=#

# Set parameters
N = 16
Δt = 10*60; # 10 minutes

# define the forward map
zᵖ = zeros(N)
start = argmin(abs.(les.t .- 86400 * 0.25)) #start at a quarter of a day
#start = 1
#subsample = start:2:length(les.t)
𝒢_n = CoreFunctionality.closure_free_convection_flexible(N, Δt, les,
                                              subsample = subsample, grid = zᵖ, power = 1.0)

NN = sqrt(les.α * les.g * les.bottom_T)
tˢ = collect(subsample) #indices for simulation time
#𝑪_n = [1e-4, 3.5 * 1.0, 10.0, 0.0, 0.375, NN]
𝑪_n = [1e-4, 3.5, 7.5, 0.0, 0.4, NN]
#𝑪_n = [9.414446198612705e-5, 3.681665772608741, 7.283972927371181, 1.1102230246251565e-16, 0.52436119436052, NN]
if optimal
    resolution_label = "_res_" * string(N)
    filename2 = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
    mcmc_data = jldopen(filename2, "r")
    chain = mcmc_data["𝑪"]
    e1 = mcmc_data["ε"]
    e2 = mcmc_data["proposal_ε"]
    close(mcmc_data)
    𝑪_n = chain[:, argmin(e1)]
    println(𝑪_n)
end



kpp_T = 𝒢_n(𝑪_n)

minT = minimum(les.T⁰)
maxT = maximum(les.T⁰)



if animation
    for j in 1:50:finish
        jˢ = tˢ[j]
        days = @sprintf("%.1f", les.t[jˢ]/86400)
        Tˢ = les.T[:,jˢ]
        Tᵖ = kpp_T[:,j]
        z  = les.z

        p1 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT))
        scatter!( Tᵖ , zᵖ, label = "KPP (h²N²)", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
        display(p1)
    end
end
jˢ = tˢ[j]
days = @sprintf("%.1f", les.t[jˢ]/86400)
Tˢ = les.T[:,jˢ]
Tᵖ = kpp_T[:,j]
z  = les.z

p2 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT))
scatter!( Tᵖ , zᵖ, label = "KPP (h²N²)", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
display(plot(p1,p2))

###
ℒᵗ_n = CoreFunctionality.closure_T_nll(𝒢_n, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )
ℒᵗ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )
ll = sqrt.(ℒᵗ(𝑪))
ll_n = sqrt.(ℒᵗ_n(𝑪_n))

plot(les.t[subsample] ./ 86400, ll, label = "old")
plot!(les.t[subsample] ./ 86400, ll_n, label = "new")


###
# Test Mixed Layer Depth Forward Map
optimal = true
# Get LES data
wd = pwd()
case = cases[1]
filename = wd * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

# Set parameters
N = 16
Δt = 10*60; # 10 minutes

# define the forward map
zᵖ = zeros(N)
start = argmin(abs.(les.t .- 86400 * 0.25)) #start at a quarter of a day
#start = 1
subsample = start:1:length(les.t)
𝒢 = CoreFunctionality.closure_free_convection_ml_depth(N, Δt, les,  subsample = subsample, grid = zᵖ)

tˢ = collect(subsample) #indices for simulation time
finish = floor(Int, length(tˢ) / 1)

# compute the forward map
𝑪 = [0.1, 4.33, 1.36, 3.19 ] #/ 3.19 * 4.32]
#𝑪 = [0.48944664115621117, 5.606733551864946, 0.06980180915616958, 11.131016010149851]
if optimal
    resolution_label = "_res_" * string(N)
    filename2 = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
    mcmc_data = jldopen(filename2, "r")
    chain = mcmc_data["𝑪"]
    e1 = mcmc_data["ε"]
    e2 = mcmc_data["proposal_ε"]
    close(mcmc_data)
    𝑪 = chain[:, argmin(e1)]
    println(𝑪)
end

h = 𝒢(𝑪)

days = @sprintf("%.1f", les.t[jˢ]/86400)

p1 = plot( les.t[subsample] ./86400 , h[1,:], legend = false , xlabel= "time [days]", ylabel = "Mixed Layer Depth [m]", title = "Mixed Layer Depth Evolution", ylims = (0,80))
display(p1)
