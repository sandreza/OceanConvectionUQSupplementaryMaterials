
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

struct RandomPlugin{S,T,V,U}
    priors::S
    fcalls::T
    seed::V
    progress::U
end

function RandomPlugin(priors, fcalls::Int; seed = 1234, progress = true)
    return RandomPlugin(priors, fcalls, 1234, true)
end



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
        return losses, vals
    end
end

## Example
# Create Prior
lower = [0.0,0.0,0.0,0.0]
upper = [1.0, 8.0, 6.0, 5.0/0.3]
priors = Uniform.(lower, upper)
# Determine number of function calls
functioncalls = 212

# Define Method
method = RandomPlugin(priors, functioncalls)
minval = optimize(ℒ, method)
