for j in 1:9
    filename = "plot" * string(j) * ".jl"
    println("running " * filename)
    include(filename)
end

include("extra_plots.jl")
include(pwd() * "/sandbox/sandbox_les.jl")
