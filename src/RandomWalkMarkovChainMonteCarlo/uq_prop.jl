"""
propagate_uncertainty(𝑪, 𝒢; field_range = [], filename = [], freq = 100)

# Description
- propagate uncertainty of forward map with respect to parameter probability distributions. Assumes only one field for the forward map

# Arguments
- `𝑪`: chain of parameters
- `𝒢`: forward map

# Keyword arguments
- `field_range`: range for partitioning the histogram
- `filename`: name of file for output
- `freq`: how often to output update
"""
function propagate_uncertainty(𝑪, 𝒢; field_range = [], filename = [], freq = 100)
    # find sizes
    Φ = 𝒢(𝑪[:,1])
    nz, nt = size(Φ)
    nℰ = length(𝑪[1,:])
    if isempty(field_range)
        ϕmin = minimum(Φ)
        ϕmax = maximum(Φ)
        Δϕ = (ϕmax - ϕmin) / 1000
        ϕrange = collect(ϕmin:Δϕ:ϕmax)
    else
        ϕrange = field_range
    end
    #first step
    h1 = []
    for i in 1:nt
        push!(h1,[])
    end
    for k in 1:nt
        ϕ = Φ[:,k]
        for j in 1:nz
            tmp = fit(Histogram, [ϕ[j]], ϕrange)
            push!(h1[k],tmp)
        end
    end

    if !isempty(filename)
        @save filename h1
    end

    # next steps
    param_length = 2:nℰ
    for i in param_length
        Φ = 𝒢(𝑪[:,i])  #forward map
        for k in 1:nt
            ϕ = Φ[:,k]
            for j in 1:nz
                tmp = fit(Histogram, [ϕ[j]], ϕrange)
                merge!(h1[k][j], tmp)
            end
        end
        if (i%freq) == 0
            if !isempty(filename)
                @save filename h1
            end
            println("done with " * string(i))
        end
    end
    if !isempty(filename)
        @save filename h1
    end
    return h1
end
