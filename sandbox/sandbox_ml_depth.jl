include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

# calculate mixed layer depth uncertainty
calc_prior = false
const number_of_ensembles = 10^6
const skip = 100
for case in cases[1:1]
    for resolution in resolutions[1:1]
        # construct filename
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        # load les
        les = CoreFunctionality.OceananigansData(filename)
        # construct default loss function
        N = resolution[1]
        Δt = resolution[2]
        # define the forward map
        zᵖ = zeros(N)
        #calculate every hour
        subsample = 1:6:length(les.t)
        # define the forward map
        𝒢 = CoreFunctionality.closure_free_convection_ml_depth(N, Δt, les, subsample = subsample, grid = zᵖ)
        println("-------------------")
        println("For case $case ")
        println("and resolution " * string(resolution[1]))
        resolution_label = "_res_" * string(resolution[1])
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
        mcmc_data = jldopen(filename, "r")
        𝑪 = mcmc_data["𝑪"]
        close(mcmc_data)
        Φ = 𝒢(𝑪)
        ϕmin = 0.0
        ϕmax = 100.0
        Δϕ = (ϕmax - ϕmin) / 1000
        ϕrange = collect(ϕmin:Δϕ:ϕmax)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_uncertainty_propagation_ml.jld2"
        CoreFunctionality.propagate_uncertainty(𝑪[:,1:skip:number_of_ensembles], 𝒢, field_range = ϕrange, filename = filename)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_domain.jld2"
        println("done with posterior")
        if calc_prior
            # Now do it for the prior distribution
            filename = pwd() * "/mcmc_data/" * "prior" * "_mcmc.jld2"
            mcmc_data = jldopen(filename, "r")
            𝑪 = mcmc_data["𝑪"]
            close(mcmc_data)
            filename = pwd() * "/mcmc_data/" * case * resolution_label * "_prior" * "_uncertainty_propagation_ml.jld2"
            CoreFunctionality.propagate_uncertainty(𝑪[:,1:skip:number_of_ensembles], 𝒢, field_range = ϕrange, filename = filename)
            println("done with prior")
        end
    end
end

###
save_figures = false
case = cases[1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)
resolution = resolutions[1]
resolution_label = "_res_" * string(resolution[1])
# get posterior data
# filename = pwd() * "/mcmc_data/" * case * resolution_label * "_prior_uncertainty_propagation_ml.jld2"
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_uncertainty_propagation_ml.jld2"
histogram_data = jldopen(filename, "r")
h1 = histogram_data["h1"]
close(histogram_data)

m = length(h1[100][1].weights)
n = length(h1)
mat = randn(m,n)
for i in 1:1:193
    ind = i
    normalization_constant = sum(h1[ind][1].weights)
    tmp = h1[ind][1].weights / normalization_constant
    @. mat[:,i] = tmp
end
ϕmin = 0.0
ϕmax = 100.0
Δϕ = (ϕmax - ϕmin) / 1000
ϕrange = collect(ϕmin:Δϕ:ϕmax)
ϕrange = (ϕrange[2:end] + ϕrange[1:end-1])./2
trange = les.t[1:6:length(les.t)] ./ 86400
p1 = heatmap(trange, ϕrange, log.(mat .+ eps(1.0))./log(10), ylabel = "Mixed Layer Depth [m]", xlabel = "days", title = "Mixed Layer Depth Uncertainty", clims = (-2.25, -1), ylims = (0, 75))
display(p1)
savefig(p1, pwd() * "/figures/ml_figure.png")
println("at time ")
println(trange[end-24])
println("days")
bools = h1[end-30][1].weights .> 8

seen = ϕrange[bools]
println(seen)
max_seen = maximum(seen[1:end])
min_seen = minimum(seen[1:end])
println("The largest mixed layer depth is $max_seen")
println("The smallest mixed layer depth is $min_seen")

timelabel = @sprintf("%.1f", trange[end])
tmp = ones(length(h1[end][1].weights))
p1 = histogram(ϕrange, weights = h1[end][1].weights, bins = 200, norm = true, xlims = (65, 75), ylabel = "probability", xlabel = "mixed layer depth", title = "Mixed layer depth uncertainty at " * timelabel * " days", legend = false)
if save_figures
    savefig(p1, pwd() * "/figures/ml_figure_alternative.png")
end
