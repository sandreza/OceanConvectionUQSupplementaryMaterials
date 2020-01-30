for j in 1:9
    filename = "plot" * string(j) * ".jl"
    println("running " * filename)
    include(filename)
end
