using Plots, Printf, Statistics, Optim, JLD2, Random

# Define cases to loop over
cases = ["high_res_general_strat_16", "high_res_general_strat_64"]
for j in 1:32
    # Get generic LES data
    case_number = string(j)
    case = "general_strat_" * case_number
    push!(cases,case)
end

mcmc_iterations = [10^5, 10^5]
for j in 1:32
    mcmc_iteration = 10^4
    push!(mcmc_iterations,mcmc_iteration)
end

cases2 = []
for j in 5:9
    # Get generic LES data
    case_number = string(j)
    case = "general_surface_forcing_" * case_number
    push!(cases2,case)
end
for j in 10:10:50
    # Get generic LES data
    case_number = string(j)
    case = "general_surface_forcing_" * case_number
    push!(cases2,case)
end

resolutions = [ (16, 10 * 60), (8, 20 * 60), (32, 5 * 60), (64, 2.5 * 60)]

function a_quantile(array, prob)
    m,n = size(array)
    quantiles = zeros(m)
    for i in 1:m
        quantiles[i] = quantile(array[i,:], prob)
    end
    return quantiles
end

# Define parameter dictonary
names = ["Surface Layer Fraction", "Nonlocal Amplitude", "Diffusivity Amplitude", "Mixing Parameter"]
values = [1,2,3,4]
parameter_dictionary = Dict(zip(values, names))

# define constants
default_𝑪 = [0.1, 6.33, 1.36, 3.19]
σ = default_𝑪 * 0.1
left_bounds = [0.0, 0.0, 0.0, 0.0]
right_bounds = [1.0, 8.0, 6.0, 16.0]


# define default closure
"""
closure_default_loss_function(filename ;N=16, Δt = 10*60, start_days = 0.25, subsample_parameter = 6)

# Description
- Returns default loss function.

# Arguments
- `filename`: file from which to load the LES data

# Keyword arguments
- `N`:(int) number of grid points in kpp
- `Δt`:(real) time step size in KPP
- `start_days`:(real) day to start simulation
- `subsample parameter`:(int) how often to subsample data file, 6 corresonds to hour
- `series`: (bool), false means only return the point-wise maximum in time

# Return
- `ℒ`: loss function. Given parameters 𝑪 outputs number corresponding to goodness of fit.
"""
function closure_default_loss_function(filename ;N=16, Δt = 10*60, start_days = 0.25, subsample_parameter = 6, series = false)
    les = CoreFunctionality.OceananigansData(filename)
    # define the forward map
    zᵖ = zeros(N)
     #start at a quarter of a day
    seconds_in_a_day = 86400
    start = argmin(abs.(les.t .- seconds_in_a_day  * start_days))
    #start = 1, 6 hours in a day
    subsample = start:subsample_parameter:length(les.t)
    # define the forward map
    𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample)
    # define the loss function
    ℒ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=series, power = 2, f1 = mean, f2 = maximum )
    return ℒ
end


# define fixed time
"""
closure_fixed_time_loss_function(filename; N=16, Δt = 10*60, final_day = 86400)

# Description
- Returns default loss function.

# Arguments
- `filename`: file from which to load the LES data

# Keyword arguments
- `N`:(int) number of grid points in kpp
- `Δt`:(real) time step size in KPP
- `start_days`:(real) day to start simulation
- `final_day`:(real). The day to calculate the error
- `series`: (bool), false means only return the point-wise maximum in time

# Return
- `ℒ`: loss function. Given parameters 𝑪 outputs number corresponding to goodness of fit.
"""
function closure_fixed_time_loss_function(filename; N=16, Δt = 10*60, final_day = 1.0, series = false)
    les = CoreFunctionality.OceananigansData(filename)
    # define the forward map
    zᵖ = zeros(N)
     #start at a quarter of a day
    seconds_in_a_day = 86400
    #start = 1, 6 hours in a day
    indmin = argmin(abs.(les.t .- final_day * seconds_in_a_day))
    subsample = indmin:indmin
    # define the forward map
    𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample)
    # define the loss function
    ℒ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=series, power = 2, f1 = mean, f2 = maximum )
    return ℒ
end


# define default closure
"""
closure_flexible_loss_function(filename ;N=16, Δt = 10*60, start_days = 0.25, subsample_parameter = 6)

# Description
- Returns flexible loss function with more parameters.

# Arguments
- `filename`: file from which to load the LES data

# Keyword arguments
- `N`:(int) number of grid points in kpp
- `Δt`:(real) time step size in KPP
- `start_days`:(real) day to start simulation
- `subsample parameter`:(int) how often to subsample data file, 6 corresonds to hour
- `series`: (bool), false means only return the point-wise maximum in time
- `power`: a keyword argument for the forward map

# Return
- `ℒ`: loss function. Given parameters 𝑪 outputs number corresponding to goodness of fit.
"""
function closure_flexible_loss_function(filename ;N=16, Δt = 10*60, start_days = 0.25, subsample_parameter = 6, series = false, power = 0.0)
    les = CoreFunctionality.OceananigansData(filename)
    # define the forward map
    zᵖ = zeros(N)
     #start at a quarter of a day
    seconds_in_a_day = 86400
    start = argmin(abs.(les.t .- seconds_in_a_day  * start_days))
    #start = 1, 6 hours in a day
    subsample = start:subsample_parameter:length(les.t)
    # define the forward map
    𝒢 = CoreFunctionality.closure_free_convection_flexible(N, Δt, les, subsample = subsample, power = power)
    # define the loss function
    ℒ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=series, power = 2, f1 = mean, f2 = maximum )
    return ℒ
end


"""
uq_prop_mat(h1; index = [])

# Description
- From histogram data for the uncertainty propagation generate a matrx

# Arguments

- `h1`: histogram

# Keyword arguments
- `index`: which index to grab data from

"""
function uq_prop_mat(h1; index = [])
    if isempty(index)
        kk = length(h1)
    else
        kk = index
    end
    h1_slice = h1[kk] # get the last one
    m = length(h1_slice[end].weights)
    n = length(h1_slice)
    mat = randn(m,n)
    normalization_constant = sum(h1[1][1].weights)
    for j in 1:n
        tmp = h1_slice[j].weights ./ normalization_constant .+ eps(1.0)
        @. mat[:,j] = tmp
    end
    return mat
end
