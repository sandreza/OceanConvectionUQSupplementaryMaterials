include("../src/LocalOceanUQSupplementaryMaterials.jl")
using Plots, Printf, Statistics

# Get LES data
wd = pwd()
filename = wd * "/LES/general_strat_16_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

# Set parameters
N = 16
Δt = 10*60; # 10 minutes

# define the forward map
zᵖ = zeros(N)
start = argmin(abs.(les.t .- 86400 * 0.25)) #start at a quarter of a day
#start = 1
subsample = start:2:length(les.t)

# define the forward map
𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les,
                                              subsample = subsample, grid = zᵖ)




# define the loss function
ℒ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=false, power = 2, f1 = mean, f2 = maximum )

# define the loss function at each moment in time
ℒᵗ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )



###
𝑪¹ = [0.1, 6.33, 3*1.36, 3.19*2]
𝑪² = [0.1, 6.33, 1.36, 3.19]
ℒ(𝑪¹)
ℒ(𝑪²)
loss1 = ℒᵗ(𝑪¹)
loss2 = ℒᵗ(𝑪²)
kpp_T1 = 𝒢(𝑪¹)
kpp_T2 = 𝒢(𝑪²)

days = les.t[subsample] ./ 86400
plot(days, loss1, label = "loss 1", legend = :bottomright)
plot!(days, loss2, label = "loss 2", xlims = "time", ylims = "celcius squared")

tˢ = collect(subsample) #indices for simulation time

minT = minimum(les.T⁰)
maxT = maximum(les.T⁰)
finish = floor(Int, length(tˢ) / 1)
for j in 1:20:finish
    jˢ = tˢ[j]
    days = @sprintf("%.1f", les.t[jˢ]/86400)
    loss_value1 = @sprintf("%.1e", loss1[j])
    loss_value2 = @sprintf("%.1e", loss2[j])

    Tˢ = les.T[:,jˢ]
    T¹ = kpp_T1[:,j]
    T² = kpp_T2[:,j]
    z  = les.z

    p1 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "loss = " * loss_value1, xlims = (minT,maxT))
    scatter!( T¹ , zᵖ, label = "KPP 1", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", xlims = (minT,maxT))

    p2 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", xlims = (minT,maxT))
    scatter!( T² , zᵖ, label = "KPP 2", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = " loss = " * loss_value2, xlims = (minT,maxT))
    #display(plot(p1))
    display(plot(p1, p2))
end
