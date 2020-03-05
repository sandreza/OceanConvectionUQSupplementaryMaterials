for j in 1:4
    filename = "plot" * string(j) * ".jl"
    println("running " * filename)
    include(filename)
end

include("extra_plots")
