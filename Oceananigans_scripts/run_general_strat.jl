const temp_array = [1]

for i in 1:1:32
  temp_array[1] = i
  println("looking at " * string(i))
  include("general_strat.jl")
end
