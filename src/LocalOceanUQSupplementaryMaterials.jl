module CoreFunctionality

using   JLD2,
        OceanTurb,
        Statistics,
        Distributions,
        Random,
        Revise,
        StatsBase

export
        torus,
        OceananigansData,
        closure_free_convection,
        closure_free_convection_flexible,
        closure_free_convection_ml_depth,
        closure_T_nll,
        closure_proposal,
        markov_chain,
        optimize,
        optimize_and_estimate_proposal,
        propagate_uncertainty


include("./RandomWalkMarkovChainMonteCarlo/rwmcmc.jl")
include("./RandomWalkMarkovChainMonteCarlo/optimize.jl")
include("./RandomWalkMarkovChainMonteCarlo/uq_prop.jl")
include("./Hook/o2o.jl")
include("./ForwardMap/fm.jl")
include("./LossFunctions/loss_functions.jl")


end
