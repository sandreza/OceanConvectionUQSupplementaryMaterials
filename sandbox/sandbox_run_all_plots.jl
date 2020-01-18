# runs all the plots

for i in 3:9
    println("generating figures in " * string(i))
    include("../figure_scripts/plot"*string(i)*".jl")
end
