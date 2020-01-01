

# Defines several functions useful for performing a random walk


"""
accept_reject(Δℒ)

# Description

- Determines the accept or reject criteria for the Monte Carlo method.

# Input: Δℒ

- `Δℒ`: (scalar) Difference of negative log likehood functions

# Output

- Boolean Value: True or False

"""
accept_reject(Δℒ) = log(rand(Uniform(0, 1))) < Δℒ

"""
markov_link(nll, 𝑪, ε, proposal)

# Description

- Takes a single step in the random walk markov chain monte carlo algorithm and outputs proposal parameters, new parameters, and the evaluate of the loss function

# Arguments

- `nll`: The negative log-likelihood function. In the absence of priors this becomes a loss function
- `𝑪`: (array), current parameter
- `ε`: (scalar), ε = nll(𝑪). The value of negative log-likelihood of the current parameter
- `proposal`: (function), determines the proposal step

# Return

- `new_𝑪`: The value of the accepted 𝑪
- `new_ε`: value of nll(new_𝑪)
- `proposal_𝑪`: The 𝑪 from the "proposal step". Was either rejected or accepted.
- `proposal_ε`: value of nll(test_𝑪)
"""
function markov_link(nll, 𝑪, ε, proposal)
    proposal_𝑪 = proposal(𝑪)
    proposal_ε = nll(proposal_𝑪)
    Δε = (ε - proposal_ε)
    if accept_reject(Δε)
        new_ε = proposal_ε
        new_𝑪 = proposal_𝑪
    else
        new_ε = ε
        new_𝑪 = 𝑪
    end
    return new_𝑪, new_ε, proposal_𝑪, proposal_ε
end



"""
markov_chain_with_save(nll, init_𝑪, proposal, nt, filename, freq)

# Description

- A random walk that computes the posterior distribution

# Arguments

- `nll`: The negative log-likelihood function. In the absence of priors this becomes a loss function
- `init_𝑪`: (Array), initial parameter values
- `proposal`: (function), proposal function for MCMC
- `nt`: (Int) number of markov chain monte carlo steps
- `perturb`: a function that performs a perturbation of 𝑪

# Keyword Arguments
- `filename`: name for output file in JLD2 format
- `freq`: how often to save output (in terms of iterations)
- `verbose`: (bool), if true then print current optimal parameters

# Return

- `param`: The matrix of accepted parameters in the random walk
- `ε`: The array of errors associated with each step in param chain

"""
function markov_chain(nll, initial_𝑪, proposal, nt;
                      filename = [], freq = 1, verbose = false)
    𝑪 = ones(length(initial_𝑪),nt+1)
    @. 𝑪[:,1] = initial_𝑪
    proposal_𝑪 = copy(𝑪)
    ε = ones(nt+1)
    proposal_ε = copy(ε)
    ε[1] = nll(initial_𝑪)
    for i in 1:nt
        new_𝑪, new_ε, proposed_𝑪, proposed_ε = markov_link(nll, 𝑪[:,i], ε[i], proposal)
        @. 𝑪[:,i+1] = new_𝑪
        ε[i+1] = new_ε
        @. proposal_𝑪[:,i+1] = proposed_𝑪
        proposal_ε[i+1] = proposed_ε
        if i%freq==0
            println("saving index " * string(i))
            if !isempty(filename)
                @save filename ε 𝑪 proposal_ε proposal_𝑪
            end
            if verbose==true
                indmin = argmin(ε[1:i])
                println("The current optimal parameters are")
                println(𝑪[:,indmin])
                println("The loss function is " * string(ε[indmin]))
                tmpstrng = string(ε[1] / ε[indmin] )
                println("This is an improvement of " * tmpstrng)
                acceptance_rate = sum(ε[1:i] .== proposal_ε[1:i]) / length(ε[1:i])
                println("The current acceptance rate is $acceptance_rate")
            end
        end
    end
    return 𝑪, ε
end

"""
torus(x, a, b)

# Description

- Takes x ∈ ℝ and outputs torus(x) ∈ [a, b] in a periodic way.
- If a particle is moving to the right then it will pop from b to the point a

# Arguments: x, a, b

- `x`: (scalar). Current location of particle
- `a`: (scalar). left endpoint of interval
- `b`: (scalar). right endpoint of interval

# Output

-  `y`: (scalar). a value in the interval [a,b]
"""
torus(x::Number, a::Number, b::Number) = (((x-a)/(b-a))%1 - 0.5 * (sign((x-a)/(b-a)) - 1) )*(b-a) + a

"""
torus(x, a, b)

# Description

- Takes x ∈ ℝⁿ and outputs torus(x) ∈ ∏[aⁿ, bⁿ] in a periodic way.
- If a particle is moving to the right then it will pop from one part of the box to the oher

# Arguments: x, a, b

- `x`: (array). Current location of particle
- `a`: (array). left endpoint of tensor product interval
- `b`: (array). right endpoint of tensor product interval

# Output

-  `y`: (array). a value in the interval ∏[aⁿ, bⁿ]
"""
function torus(x::AbstractArray, a::AbstractArray, b::AbstractArray)
    N = length(x)
    y = zeros(N)
    for i in 1:N
        y[i] = torus(x[i], a[i], b[i])
    end
    return y
end


"""
closure_proprosal(covariance = Σ; left_bounds = [], right_bounds = []))

# Description

- Constructs a proposal for the Monte Carlo method.

# Arguments

- `covariance`: (vector) proposal parameter

# Keyword Arguments

- `left_bounds`: (array), left bounds for parameters
- `right_bounds`: (array), right bounds for parameters

# Output:

- `proposal`: (function), a function that outputs the proposal parameter

"""
function closure_proposal(Σ; left_bounds = [], right_bounds = [])
    perturbation = MvNormal(Σ)
    function proposal(𝑪)
        proposal_𝑪 = copy(𝑪)
        proposal_𝑪 .+= rand(perturbation)
        # limit ranges for the parameters
        if isempty(left_bounds)
            return proposal_𝑪
        else
            return torus(proposal_𝑪, left_bounds, right_bounds)
        end
        return nothing
    end
    return proposal
end
