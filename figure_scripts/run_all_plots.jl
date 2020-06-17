for j in 1:9
    filename = "plot" * string(j) * ".jl"
    println("running " * filename)
    include(filename)
end

include("extra_plots.jl")

# This has the mcmc example
include("mcmc_examples.jl")

# This has the toy mcmc examples
include(pwd() * "/sandbox/sandbox_les.jl")
