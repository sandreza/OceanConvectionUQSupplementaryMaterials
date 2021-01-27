include("../src/LocalOceanUQSupplementaryMaterials.jl")


using Plots, Printf, Statistics, LinearAlgebra, OceanTurb
using Random, Distributions
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

nt = length(subsample)
coarseT = zeros(N, nt)
for i in 1:nt
    coarseT[:,i] .= CoreFunctionality.avg( les.T[:, subsample[i]], N)
end
# define the loss function
ℒ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=false, power = 2, f1 = mean, f2 = maximum )

# h 
params = [0.1, 6.33, 1.36, 3.19]
ℒ(params)

function gradient(L, params; δ = 1e-2 .* params)
    ∇L = zeros(length(params))
    e = I + zeros(length(params), length(params))
    for i in eachindex(params)
        up = L(params + δ[i] * e[i,:])
        down = L(params - δ[i] * e[i,:])
        ∇L[i] = (up - down)/(2*δ[i])
    end
    return ∇L  
end

function gradient(L, params; δ = 1e-4 .* params)
    ∇L = zeros(length(params))
    e = I + zeros(length(params), length(params))
    Lcurrent = L(params)
    for i in eachindex(params)
        up = L(params + δ[i] * e[i,:])
        ∇L[i] = (up - Lcurrent)/(δ[i])
    end
    return ∇L  
end

params = [0.1,0.1,0.1,0.1]
∇L = gradient(ℒ, params; δ = 1e-4 .* params)
L0 = ℒ(params)
α1 = 1e-3 / L0
L1 = ℒ(params - α1 * ∇L)
α2 = 2e-3 / L0
L2 = ℒ(params - α2 * ∇L)
α3 = 3e-3 / L0
L3 = ℒ(params - α3 * ∇L)

αm1 = -1e-3 / L0
LM1 =  ℒ(params - αm1 * ∇L)

αm2 = -2e-3 / L0
LM2 =  ℒ(params - αm2 * ∇L)


α4 = 4e-3 / L0
L4 = ℒ(params - α4 * ∇L)


αs = range(-0e-3 / L0, 1e-3/L0, length = 10)
LS = [ℒ(params - α * ∇L) for α in αs]
plot(αs, LS)
plot([LM2, LM1, L0, L1, L2, L3, L4])
##
lower = [0.0,0.0,0.0,0.0]
upper = [1.0, 8.0, 6.0, 5.0/0.3]
fcalls = 212
Random.seed!(1234)
guessparams = [0.1, 0.1, 0.1, 0.1]
bestloss = [ℒ(guessparams)]
bestparams = copy(guessparams)
for i in 1:fcalls
    guessparams = rand.(Uniform.(lower, upper))
    currentloss = ℒ(guessparams)
    if currentloss < bestloss[1]
        bestloss[1] = currentloss
        bestparams .= guessparams
        println("found something better at ", i)
    end
end
println("best after guessing")
println((bestparams, bestloss))
linesearches = 10
αs = collect(range(-0e-3 / bestloss[1], 1e-3 / bestloss[1], length = 10))
for i in 1:linesearches
    ∇L = gradient(ℒ, bestparams; δ = 1e-4 .* upper)  
    LS = [ℒ(bestparams - α * ∇L) for α in αs]
    b = argmin(LS)
    println("best after line search ", i)
    bestparams .= bestparams - αs[b] * ∇L
    if αs[b] == 0
        αs .= 0.5 .* αs
        println("staying the same at ")
    end
    println((bestparams, LS[b]))
end