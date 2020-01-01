using MCMCDiagnostics
using JLD2
include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

case = cases[1]
resolution = resolutions[1]
N = resolution[1]
resolution_label = "_res_" * string(resolution[1])
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
resolution_label = "_res_" * string(N)
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
close(mcmc_data)

a_ratio = sum(e1 .== e2) / length(e1)
println("acceptance ratio")
println(a_ratio)
# check effecitve sample size

ess = randn(4)
for i in 1:4
    # one million is a bit much
    x1 = chain[i,1:100000]
    variance_x1 = var(x1)
    ess[i] = effective_sample_size(x1, variance_x1)
end
println("effectiive sample size")
println(ess)
