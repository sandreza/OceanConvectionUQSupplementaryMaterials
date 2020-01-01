include("../src/LocalOceanUQSupplementaryMaterials.jl")
using Plots, Printf

wd = pwd()
filename = wd * "/LES/general_strat_1_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

plot(les.T⁰, les.z, xlabel = "Temperature [C]", ylabel = "Depth [m]")
plot!(legend = false)



minT = minimum(les.T⁰)
maxT = maximum(les.T⁰)

start  = argmin(abs.(les.t .- 86400 * 0.25))
finish = floor(Int, length(les.t) / 1)
for j in start:5:finish
    days = @sprintf("%.1f", les.t[j]/86400)

    Tˢ = les.T[:,j]
    z  = les.z

    p1 = plot( Tˢ , z, label = "LES", legend = :bottomright , xlabel= "Temperature [C]", ylabel = "Depth [m]", title = "Vertical Profile at t = " * days * " days", xlims = (minT,maxT))
    display(p1)
end
