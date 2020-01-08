using LaTeXStrings

const celsius = L"[$^\circ$C]";
itime = L"[s^{-2}]"
# utils for plots
"""
get_chain(case, N)

# Description
- Get markov chain data from casenumber and gridpoint label

# Arguments
- 'case': (string). case identifier
- 'N': (int). number of gridpoints used in mcmc

# Return
- `chain`: (array). markov chain data
- `e1`: (vector). loss function values associated with markov chain data
- `e2`: (vector). loss function values of the proposal
"""
function get_chain(case, N)
        resolution_label = "_res_" * string(N)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
        mcmc_data = jldopen(filename, "r")
        chain = mcmc_data["𝑪"]
        e1 = mcmc_data["ε"]
        e2 = mcmc_data["proposal_ε"]
        close(mcmc_data)
        return chain, e1, e2
end

"""
get_optimal_guess(case, N)

# Description
- Get optimization data from casenumber and gridpoint label

# Arguments
- 'case': (string). case identifier
- 'N': (int). number of gridpoints used in mcmc

# Return
- `parameter`: (array). optimal parameter data
- `loss`: (vector). loss function values associated with data
"""
function get_optimal_guess(case, N)
        resolution_label = "_res_" * string(N)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_optima.jld2"
        opt_data = jldopen(filename, "r")
        param = opt_data["parameter"]
        loss = opt_data["loss"]
        close(opt_data)
        return param, loss
end


"""
calculate_partial_statistics(chain)

# Description
- Calculates the mean for subsets of the data, e.g., μ[1,n] takes into account the first n+1 data points in the markov chain for the first state variable but neglects the rest

# Arguments
- `chain`: (array), markove chain

# Return
- `μ`: (array). partial means of every state
- `σ`: (array). partial standard deviations of every state
"""
function calculate_partial_statistics(chain)
    nt = length(chain[1, :])
    μ = zeros(4, nt)
    σ = zeros(4, nt)

    @. μ[:, 1] = chain[:,1]

    for j in 2:nt
            for s in 1:4
                    δμ = (chain[s,j] - μ[s, j-1]) / j
                    μ[s, j] =  μ[s, j-1] + δμ
                    δσ =  ((chain[s,j] - μ[s, j])^2 - σ[s, j-1]) / (j-1)
                    σ[s, j] =  σ[s, j-1] + δσ + δμ^2
            end
    end
    return μ, sqrt.(σ)
end

function calculate_partial_statistics2(chain)
    nt = length(chain[1, :])
    μ = zeros(4, nt-1)
    σ = zeros(4, nt-1)

    for j in 1:(nt-1)
            for s in 1:4
                    μ[s, j] = mean(chain[s,1:(j+1)])
                    σ[s, j] = std(chain[s,1:(j+1)])
            end
    end
    return μ, σ
end

"""
plot_partial_statistics(μ, σ; ind = 1)

# Description

- Plots partial statistcs calculated from function calculate_partial_statistics()

# Arguments
- `μ`: (array). partial means of every state
- `σ`: (array). partial standard deviations of every state

# Keyword Arguments
- `index`: (int). The state that we would like to plot

"""
function plot_partial_statistics(μ, σ; ind = 1)
        pμ = plot(μ[ind,:], legend = false, xlabel = "iteration", ylabel = "partial mean", title = parameter_dictionary[ind], grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
        pσ = plot(σ[ind,:], legend = false, xlabel = "iteration", ylabel = "partial std", title = parameter_dictionary[ind], grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
        return pμ, pσ
end

"""
marginal_pdfs(chain, left_bounds, right_bounds, parameter_dictionary)

# Description
- Generate marginal distribution pdfs for each state variable

# Arguments
- `chain`: (array). mcmc result
- `left_bounds`: (vector), lower bound for plot
- `right_bounds`: (vector), upper bound for plot
-  `parameter_dictionary`: (vector), vector of strings for labeling the state name

# Keyword Arguments

- `bins`: (int): number of bins used to construct histogram

# Return
- `p`: plot container. p[1] has the marginal distribution of the first state variable. to plot all the states at once use plot(p...)
"""
function marginal_pdfs(chain, left_bounds, right_bounds, parameter_dictionary; bins = 50)
        p = []
        m,n = size(chain)
        for i in 1:m
                index = i
                p1 = histogram(chain[index, :], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = false, normalize = true, ylabel = "pdf", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
                push!(p, p1)
        end
        return p
end

"""
marginal_pdfs(chain1, chain2, left_bounds, right_bounds, parameter_dictionary)

# Description
- Generate marginal distribution pdfs for each state variable

# Arguments
- `chain1`: (array). mcmc result of case 1
- `chain2`: (array). mcmc result of case 2
- `left_bounds`: (vector), lower bound for plot
- `right_bounds`: (vector), upper bound for plot
-  `parameter_dictionary`: (vector), vector of strings for labeling the state name

# Keyword Arguments

- `bins`: (int): number of bins used to construct histogram

# Return
- `p`: plot container. p[1] has the marginal distribution of the first state variable. to plot all the states at once use plot(p...)
"""
function marginal_pdfs(chain1, chain2, left_bounds, right_bounds, parameter_dictionary; bins = 50)
        p = []
        for i in 1:4
                index = i
                p1 = histogram(chain1[index, :], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = true, normalize = true, ylabel = "pdf", alpha = 0.4, label = "Case 1", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
                p1 = histogram!(chain2[index, :], xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = bins, legend = true, normalize = true, ylabel = "pdf", alpha = 0.4, label = "Case 2", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
                push!(p, p1)
        end
        return p
end

"""
joint_pdfs(chain, left_bounds, right_bounds, parameter_dictionary; indpairs)

# Description
- Generate joint distribution pdfs for each combination of state variables

# Arguments
- `chain`: (array). mcmc result
- `left_bounds`: (vector), lower bound for plot
- `right_bounds`: (vector), upper bound for plot
-  `parameter_dictionary`: (vector), vector of strings for labeling the state name

# Keyword Arguments
- `indpairs`: index pairs to generate pdfs, by default assumes 4
- `bins`: (int), number of bins for constructing the histogram

# Return
- `p`: plot container. p[1] has the joint distribution of two state variables. to plot all combinations at once use plot(p...)
"""
function joint_pdfs(chain, left_bounds, right_bounds, parameter_dictionary; indpairs = [[1,2], [1,3], [1,4], [2,3], [2,4], [3,4]], bins = 50)
        p = []
        for pair in indpairs
                index1 = pair[1]
                index2 = pair[2]
                p1 = histogram2d(chain[index1, :], chain[index2, :], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
                push!(p, p1)
        end
        return p
end


"""
combine(chain1, chain2)

# Description
- Combine two pdfs together.

# Arguments
- `chain1`: (array). First markov chain
- `chain2`: (array). Second markov chain

# Return
- `new_chain`: (array). new markov chain
"""
function combine(chain1, chain2)
        ns, ne1 = size(chain1)
        ns, ne2 = size(chain2)
        new_chain = zeros(ns, ne1 + ne2)
        @. new_chain[:,1:ne1] = chain1
        @. new_chain[:,(ne1+1):end] = chain2
        return new_chain
end
