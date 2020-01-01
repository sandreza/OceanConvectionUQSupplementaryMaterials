include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
using Plots, Printf, Statistics, JLD2

save_figures = true

case = cases[1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

 #start at a quarter of a day
seconds_in_a_day = 86400
start = 1
#start = 1, 6 hours in a day
subsample_parameter = 1
subsample = start:subsample_parameter:length(les.t)


# forward map parameters
p = []
labels = ["Resolution 1", "Resolution 2", "Resolution 3"]
Nlist = [16, 64, 512]
Δt = [10*60, 3*60, 30]
plot()
for j in zip(Nlist, Δt, labels)
    N = j[1]
    Δt = j[2] #seconds
    zᵖ = zeros(N)
    # define the forward map
    𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample)
    # define the loss function
    ℒᵗ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )

    t = les.t ./ seconds_in_a_day
    inds = 30:length(les.t)
    loss = ℒᵗ(default_𝑪)
    p1 = plot!(t[inds], sqrt.(loss[inds]), label = j[3], legend = :bottomright, xlabel = "days", ylabel = "Error " * "[C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box )
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_1.png")
    end
    display(p1)
    push!(p,p1)
end
