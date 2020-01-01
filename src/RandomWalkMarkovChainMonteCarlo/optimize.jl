
"""
optimize(initial_𝑪, nll; nt = 1000, restart = 0, proposal = [], scale = 1, filename = [], rescale = true, freq = 1001)

# Description
- A generic optimizer using RWMCMC. It is generally better to use Optim

# Arguments

-  `initial_𝑪`:(vector) initial parameter
- `nll`:(function) negative log-likelihood. The function to minimize

# Keyword Arguments
- `nt`: (int), how many steps of the random walk to take
- `restart`: (int), restart at the optimal value this many times
- `proposal`: (function), proposal function for performing random walk
- `scale`: (real), scale for constructing default proposal
- `filename`: (string), a place to save the JLD2 file for optimization
- `rescale`: (boolean), allows one to rescale the loss function over iterations
- `freq`: how often to output progress, make this larger than nt for no output
- `verbose`: (boolean), outputs optimal values with frequence = freq

# Comments
- This is a "prep step" in mcmc

"""
function optimize(initial_𝑪, nll; nt = 10000, restart = 0, proposal = [], scale = 0.2, filename = [], rescale = true, freq = 10001, verbose = true)
    if proposal == []
        perturbation = closure_proposal(initial_𝑪 * scale)
    else
        perturbation = proposal
    end
    if rescale == true
        scale = nll(initial_𝑪)
    else
        scale = 1.0
    end
    ℒ(𝑪) = nll(𝑪) / scale
    # perform random walk
    tmp_𝑪 = copy(initial_𝑪)
    for i in 1:(restart+1)
        new_𝑪, new_ε = markov_chain(ℒ, tmp_𝑪, perturbation, nt; freq = freq, filename = filename, verbose = verbose)
        # pick out new optimal value
        optimal_index = argmin(new_ε)
        opt_𝑪 = new_𝑪[:, optimal_index]
        tmp_𝑪 = opt_𝑪
        if rescale == true
            ℒ(𝑪) = nll(𝑪) / nll(tmp_𝑪)
        end
    end
    return tmp_𝑪
end


"""
optimize_and_estimate_proposal(initial_𝑪, nll, left_bounds, right_bounds; nt = 1000, restart = 0, proposal = [], scale = 1, filename = [], rescale = true, freq = 1001)

# Description
- A generic optimizer using RWMCMC. It also tries to estimate a new proposal

# Arguments

-  `initial_𝑪`:(vector) initial parameter
- `nll`:(function) negative log-likelihood. The function to minimize
- `left_bounds`: bounds for the proposal
- `right_bounds`: bounds for the proposal

# Keyword Arguments
- `nt`: (int), how many steps of the random walk to take
- `restart`: (int), restart at the optimal value this many times
- `proposal`: (function), proposal function for performing random walk
- `scale`: (real), scale for constructing default proposal
- `filename`: (string), a place to save the JLD2 file for optimization
- `rescale`: (boolean), allows one to rescale the loss function over iterations
- `freq`: how often to output progress, make this larger than nt for no output
- `verbose`: (boolean), outputs optimal values with frequence = freq

# Comments
- This is a "prep step" in mcmc

"""
function optimize_and_estimate_proposal(initial_𝑪, nll, left_bounds, right_bounds; nt = 10000, restart = 0, proposal = [], scale = 0.2, filename = [], rescale = true, freq = 10001, verbose = true)
    if proposal == []
        perturbation = closure_proposal(initial_𝑪 * scale, left_bounds = left_bounds, right_bounds = right_bounds)
    else
        perturbation = proposal
    end
    if rescale == true
        scale = nll(initial_𝑪)
    else
        scale = 1.0
    end
    ℒ(𝑪) = nll(𝑪) / scale
    # perform random walk
    tmp_𝑪 = copy(initial_𝑪)
    Σ = randn(length(initial_𝑪),length(initial_𝑪))
    for i in 1:(restart+1)
        new_𝑪, new_ε = markov_chain(ℒ, tmp_𝑪, perturbation, nt; freq = freq, filename = filename, verbose = verbose)
        # pick out new optimal value
        optimal_index = argmin(new_ε)
        opt_𝑪 = new_𝑪[:, optimal_index]
        tmp_𝑪 = opt_𝑪
        tmp_Σ = cov(new_𝑪')
        println(Σ)
        @. Σ = tmp_Σ
        perturbation = closure_proposal(Σ, left_bounds = left_bounds, right_bounds = right_bounds)
        if rescale == true
            ℒ(𝑪) = nll(𝑪) / nll(tmp_𝑪)
        end
    end
    return tmp_𝑪, Σ
end
