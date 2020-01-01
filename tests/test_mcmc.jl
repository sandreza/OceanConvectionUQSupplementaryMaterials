include("../src/LocalOceanUQSupplementaryMaterials.jl")
using Plots, Printf, JLD2, Revise, Statistics, Random

# Define negative log-likelihood
σ = 1.0
nll(𝑪) = 1/(2σ^2) * 𝑪[1]^2
# initial parameter for MCMC
initial_𝑪 = [1.0]
# construct proposal matrix, proposal step can be undrelated to distribution
proposal = CoreFunctionality.closure_proposal([2.0])
proposal([1.0])
# determine number of steps in chain
nt = 10000
sd = pwd() * "/mcmc_data/test.jld2"
# for reproducibility
Random.seed!(1234)
CoreFunctionality.markov_chain(nll, initial_𝑪, proposal,
                               nt,  freq = 10, filename = sd)

mcmc_data = jldopen(sd, "r")
chain = mcmc_data["𝑪"]
close(mcmc_data)

C = collect(-4σ:0.1:4σ)

for j in 1:1000:length(chain[1,:])
        histogram(chain[1,1:j], normalize = true, alpha = 0.5,
                  bins = 100, label = "histogram")
        p1= plot!(C, exp.(-nll.(C))./sqrt(2 * π * σ^2),
                color = :black, linewidth = 2.0, xlims = (-4σ,4σ), ylims = (0, 1.1/sqrt(2 * π * σ^2)), xlabel = "parameter", ylabel = "pdf", label = "exact solution")
        display(p1)
end



###
# show ergodic chain
plot(chain[1,:], xlabel = "iteration", ylabel = "parameter")

# show convergence of mean and standard deviation
μ = zeros(nt-1)
standard_deviation = zeros(nt-1)
for j in 1:(nt-1)
        μ[j] = mean(chain[1,1:(j+1)])
        standard_deviation[j] = std(chain[1,1:(j+1)])
end
its = collect(1:(nt-1))
plot(log.(its), log.(abs.(μ)), xlabel = "iteration", label = "log mean error")
plot!(log.(its), log.(abs.(standard_deviation .- σ)), xlabel = "iteration", label = "log standard deviation error", legend = :bottomleft)
println("final mean is $(μ[end])")
plot!(log.(its), - 1/2 * log.(its) .+ 1.0,label = "theoretical convergence rate", width = 2)
println("final standard deviation is $(standard_deviation[end])")


###
# test non-guassian distributions
# need to keep the speration wide so that the approximation to the
# normalization constant is accurate
nll(𝑪) = 1/(2σ^2) * (𝑪[1]^2 - 4)^2
# initial parameter for MCMC
initial_𝑪 = [1.0]
# construct proposal matrix, proposal step can be undrelated to distribution
proposal = CoreFunctionality.closure_proposal([2.0])
proposal([1.0])
# determine number of steps in chain
nt = 10000
sd = pwd() * "/mcmc_data/test2.jld2"
# for reproducing data
Random.seed!(1234)
CoreFunctionality.markov_chain(nll, initial_𝑪, proposal,
                               nt,  freq = 10, filename = sd)

mcmc_data = jldopen(sd, "r")
chain = mcmc_data["𝑪"]
close(mcmc_data)

C = collect(-4σ:0.05:4σ)

for j in 1:100:length(chain[1,:])
        histogram(chain[1,1:j], normalize = true, alpha = 0.5,
                  bins = 100, label = "histogram")
        p1= plot!(C, 2.0 * exp.(-nll.(C))./sqrt(2 * π * σ^2),
                color = :black, linewidth = 2.0, xlims = (-4σ,4σ), ylims = (0, 2.1/sqrt(2 * π * σ^2)), xlabel = "parameter", ylabel = "pdf", label = "approximate solution", legend = :top)
        display(p1)
end


plot(chain[1,:], xlabel = "iteration", ylabel = "parameter")

###
# Test mulitivariate distribution
# Define negative log-likelihood

const σ = 1.0
nll(𝑪) = 1/(2σ^2) * (𝑪[1]^2 + (𝑪[2]-1)^2 + 𝑪[1]*𝑪[2])
const 𝒁1 = 2 * exp(1/(6*σ^2)) * sqrt(2π/3) * sqrt(1/σ^2)
nll_marginal1(x) = (exp(-((x * (4 + 3*x) ) / (8*σ^2)))) / 𝒁1 ;
const 𝒁2 = 2 * exp(1/(6*σ^2)) * sqrt(2π/3) * sqrt(1/σ^2)
nll_marginal2(y) = exp(-(((-2 + y)*(-2 + 3*y))/(8* σ^2))) / 𝒁2 ;
# initial parameter for MCMC
initial_𝑪 = [1.0, -1.0]
# construct proposal matrix, proposal step can be undrelated to distribution
proposal = CoreFunctionality.closure_proposal([2.0, 2.0])
# determine number of steps in chain
nt = 10000
sd = pwd() * "/mcmc_data/test3.jld2"
# for reproducibility
Random.seed!(1234)
CoreFunctionality.markov_chain(nll, initial_𝑪, proposal,
                               nt,  freq = 10, filename = sd)

mcmc_data = jldopen(sd, "r")
chain = mcmc_data["𝑪"]
close(mcmc_data)

C = collect(-4σ:0.1:4σ)

for j in 1:1000:length(chain[1,:])
        histogram(chain[1,1:j], normalize = true, alpha = 0.5,
                  bins = 100, label = "histogram")
        p1= plot!(C .- 0.5, nll_marginal1.(C .- 0.5),
                color = :black, linewidth = 2.0, xlims = (-4σ .- 0.5,4σ .- 0.5), ylims = (0, 1.1/sqrt(2 * π * σ^2)), xlabel = "parameter", ylabel = "pdf", label = "exact solution")

        histogram(chain[2,1:j], normalize = true, alpha = 0.5,
                  bins = 100, label = "histogram")
        p2 = plot!(C .+ 1.5, nll_marginal2.(C.+ 1.5),
                color = :black, linewidth = 2.0, xlims = (-4σ + 1.5,4σ + 1.5), ylims = (0, 1.1/sqrt(2 * π * σ^2)), xlabel = "parameter", ylabel = "pdf", label = "exact solution")
        display(plot(p1,p2))
end
