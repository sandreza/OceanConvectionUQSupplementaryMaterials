include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")

const number_of_ensembles = 10000
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
        𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample, grid = zᵖ)
        println("-------------------")
        println("For case $case ")
        println("and resolution " * string(resolution[1]))
        resolution_label = "_res_" * string(resolution[1])
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
        mcmc_data = jldopen(filename, "r")
        𝑪 = mcmc_data["𝑪"]
        close(mcmc_data)
        Φ = 𝒢(𝑪)
        ϕmin = minimum(Φ)
        ϕmax = maximum(Φ)
        Δϕ = (ϕmax - ϕmin) / 1000
        ϕrange = collect(ϕmin:Δϕ:ϕmax)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_uncertainty_propagation.jld2"
        CoreFunctionality.propagate_uncertainty(𝑪[:,1:skip:number_of_ensembles], 𝒢, field_range = ϕrange, filename = filename)
        # to construct grid
        𝒢(𝑪[:,1]);
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_domain.jld2"
        @save filename zᵖ ϕrange
        println("done with posterior")
        # Now do it for the prior distribution
        filename = pwd() * "/mcmc_data/" * "prior" * "_mcmc.jld2"
        mcmc_data = jldopen(filename, "r")
        𝑪 = mcmc_data["𝑪"]
        close(mcmc_data)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_prior" * "_uncertainty_propagation.jld2"
        CoreFunctionality.propagate_uncertainty(𝑪[:,1:skip:number_of_ensembles], 𝒢, field_range = ϕrange, filename = filename)
        println("done with prior")
    end
end
