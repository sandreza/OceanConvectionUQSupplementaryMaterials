include("../src/LocalOceanUQSupplementaryMaterials.jl")
using Statistics, Distributions, Random, BenchmarkTools, Optim, JLD2

# Define negative log-likelihood
const σ = 1.0
nll(𝑪) = 1/(2σ^2) * ((𝑪[1]-𝑪[2])^2 + 𝑪[1]^2)
# initial parameter for MCMC
initial_𝑪 = [1.0, 2.0]
# construct proposal matrix, proposal step can be unrelated to distribution
proposal = CoreFunctionality.closure_proposal([2.0, 2.0])
# determine number of steps in chain
nt = 100
restart = 10
rescale = true
Random.seed!(1234)
optima = CoreFunctionality.optimize(initial_𝑪, nll; nt = nt, restart = restart, proposal = proposal, scale = 1, filename = [], rescale = rescale)

println("The value of the discovered optima is $optima")


###
# test in simple case whether or not use use built-in optimization toolkit
@btime CoreFunctionality.optimize(initial_𝑪, nll; nt = nt, restart = restart, proposal = proposal, scale = 1, filename = [], rescale = rescale);

@btime optimize(nll, initial_𝑪);

res = optimize(nll,  initial_𝑪)
optim_min = Optim.minimizer(res)
println("The minima found by optim is $(optim_min)")
