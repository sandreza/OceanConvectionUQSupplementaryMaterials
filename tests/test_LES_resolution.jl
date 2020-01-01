include("../src/LocalOceanUQSupplementaryMaterials.jl")
using Plots, Printf


# Load 2 LES profiles,

wd = pwd()
filename = wd * "/LES/general_strat_16_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

wd = pwd()
filename = wd * "/LES/high_res_general_strat_16_profiles.jld2"
high_res_les = CoreFunctionality.OceananigansData(filename)


minT = minimum(les.T⁰)
maxT = maximum(les.T⁰)

start  = argmin(abs.(les.t .- 86400 * 0.25))
finish = floor(Int, length(les.t) / 1)

# plot low res and high res
for j in start:5:finish
    days = @sprintf("%.1f", les.t[j]/86400)

    Tˢ = les.T[:,j]
    hTˢ = high_res_les.T[:,j]
    z  = les.z
    hz = high_res_les.z

    p1 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT))
    p1 = plot!( hTˢ , hz, label = "High Res LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT))
    display(p1)
end


# plot difference

for j in start:5:finish
    days = @sprintf("%.1f", les.t[j]/86400)

    Tˢ = les.T[:,j]
    hTˢ = high_res_les.T[:,j]
    ΔT = CoreFunctionality.avg(hTˢ, length(Tˢ)) .- Tˢ
    z  = les.z

    p1 = plot( ΔT, z, label = "LES", legend = :bottomright , xlabel= "Systematic Bias [C]", ylabel = "Depth [m]", title = "Vertical Profile Difference at t = " * days * " days", xlims = (-0.065, 0.05))
    display(p1)
end


# another way of visualizing systematic error

Tˢ = les.T[:,end]
hTˢ = high_res_les.T[:,end]
graining_index = 0:1:7
array_maxΔT = collect(graining_index) * 1.0
array_meanΔT = collect(graining_index) * 1.0
for i in graining_index
    grain = 2^i
    println("Coarse grained bias for $grain gridpoints")
    coarse_hTˢ = CoreFunctionality.avg(hTˢ, grain)
    coarse_Tˢ = CoreFunctionality.avg(Tˢ, grain)
    coarse_z = CoreFunctionality.avg(high_res_les.z, grain)
    ΔT = coarse_hTˢ .- coarse_Tˢ
    maxΔT = maximum( abs.(ΔT) )
    meanΔT = mean( abs.(ΔT) )
    array_maxΔT[i+1] = maxΔT
    array_meanΔT[i+1] = meanΔT
    println("max is $maxΔT and the mean is $meanΔT")
    if i == 0
        p1 = scatter(ΔT, coarse_z, label = "$grain", legend = :topleft)
    else
        p1 = scatter!(ΔT, coarse_z, label = "$grain")
    end
    display(p1)
end

indices = collect(graining_index)
gridpoints = exp.(log(2) .* indices)
scatter(gridpoints, array_maxΔT, xlabel = "power -1", ylabel = "maximum bias", legend = false)
scatter!(gridpoints, array_meanΔT, xlabel = "power -1", ylabel = "mean bias", legend = false)
