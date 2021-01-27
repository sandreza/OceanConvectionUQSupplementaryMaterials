using Random, Distributions, LinearAlgebra
# numerical gradient
function gradient(L, params; δ = 1e-4 .* ones(length(params)))
    ∇L = zeros(length(params))
    e = I + zeros(length(params), length(params))
    Lcurrent = L(params)
    for i in eachindex(params)
        up = L(params + δ[i] * e[i,:])
        ∇L[i] = (up - Lcurrent)/(δ[i])
    end
    return ∇L  
end

# Define Method structs
struct RandomPlugin{S,T,V,U}
    priors::S
    fcalls::T
    seed::V
    progress::U
end

function RandomPlugin(priors, fcalls::Int; seed = 1234, progress = true)
    return RandomPlugin(priors, fcalls, 1234, true)
end

struct RandomLineSearch{I, T, B}
    linesearches::I
    linelength::I
    linebounds::T
    progress::B
    seed::I
end

function RandomLineSearch(linesearches, linelength)
    return RandomLineSearch(linesearches, linelength, (-0.1,1), true, 1234)
end
function RandomLineSearch(; linesearches = 10, linelength = 10, linebounds = (-0.1,1), progress = true, seed = 1234)
    return RandomLineSearch(linesearches, linelength, linebounds, progress, seed)
end

# Define Helper functions

# RandomPlugin
function priorloss(ℒ, fcalls, priors; seed = 1234, progress = true)
    Random.seed!(seed)
    losses = []
    vals = []
    for i in 1:fcalls
        guessparams = rand.(priors)
        push!(vals, guessparams)
        currentloss = ℒ(guessparams)
        push!(losses, currentloss)
        if progress
            println("iteration " * string(i))
        end
    end
    return losses, vals
end

# LineSearch
function randomlinesearch(ℒ, ∇ℒ, bparams; linesearches = 10, linelength = 10, linebounds = (-0.1, 1.0), progress = true, seed = 1234)
    params = copy(bparams)
    Random.seed!(seed)
    for i in 1:linesearches
        αs = [rand(Uniform(linebounds...), linelength-1)..., 0.0]
        ∇L = ∇ℒ(params) 
        LS = [ℒ(params - α * ∇L) for α in αs]
        b = argmin(LS)
        if progress
            println("best after line search ", i)
        end
        params .= params - αs[b] * ∇L
        if αs[b] == 0
            αs .= 0.5 .* αs
            if progress
                println("staying the same at ")
            end
        end
        if progress
            println((params, LS[b]))
        end
    end
    return params
end

# Define optimize functions
function optimize(ℒ, method::RandomPlugin; history = false, printresult = true)
    losses, vals = priorloss(ℒ, method.fcalls, 
                       method.priors, seed = method.seed, 
                       progress = method.progress)
    indmin = argmin(losses)
    if printresult
        println("The minimum loss is ", losses[indmin])
        println("The minimum argument is ", vals[indmin])
        println("This occured at iteration ", indmin)
    end
    if history == false
        return vals[indmin]
    else
        return vals, losses
    end
end
optimize(ℒ, initialguess, method::RandomPlugin; history = false, printresult = true) = optimize(ℒ, method::RandomPlugin; history = false, printresult = true)

function optimize(ℒ, ∇ℒ, params, method::RandomLineSearch; history = false, printresult = true)
    bestparams = randomlinesearch(ℒ, ∇ℒ, params; linesearches = method.linesearches,
                    linelength = method.linelength,
                    linebounds = method.linebounds, 
                    progress = method.progress, 
                    seed = method.seed)
    return bestparams
end

## Example
# loss function
loss(x) = (x[1] - 1)^2 + (x[2] - 2)^2 + (x[3] - 3)^2 + (x[4]-4)^2

# First construct global search
# Create Prior
lower = [0.0, 0.0, 0.0,  0.0]
upper = [2.0, 3.0, 4.0,  5.0]
priors = Uniform.(lower, upper)
# Determine number of function calls
functioncalls = 1000
# Define Method
method = RandomPlugin(priors, functioncalls)
# Optimize
minparam = optimize(loss, method)

# Next do gradient descent
# construct numerical gradient
∇loss(params) = gradient(loss, params)
# optimize choosing minimum from the global search for refinement
params = minparam
method  = RandomLineSearch(linebounds = (0, 1e-0/norm(∇loss(params))), linesearches = 20)
bestparam = optimize(loss, ∇loss, params, method)