# choose case
case = cases2[1]

# get LES
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

# define loss functions and forward maps
subsample = 1:1:length(les.t)
N = 16
Δt = 10 * 60 #seconds
zᵖ = zeros(N)
# define the forward map
𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample, grid = zᵖ)

les.top_T

les.bottom_T

-(mean(les.T,dims=1)[end] - mean(les.T,dims=1)[1] ) / (les.t[end] - les.t[1])
-les.top_T / les.L
Tᵖ = 𝒢(default_𝑪)
plot(les.t, mean(les.T,dims=1)[:], label = "LES")
plot!(les.t, mean(Tᵖ,dims=1)[:], label = "KPP")
plot!(les.t,mean(les.T,dims=1)[1] .- les.top_T .* les.t ./ les.L, label = "analytic")

exact_slope = 50 / (les.ρ * les.cᵖ) / les.L

error = abs( (mean(les.T[:,1],dims=1)[end] - mean(les.T[:,end],dims=1)[1] ) / (les.t[1] - les.t[end]) +   exact_slope) / exact_slope * 100

println("error is $error percent")



###
resolution = resolutions[3]
for i in 3:34

case = cases[i]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)
stratification = les.α * les.g * les.bottom_T
chain, tmp1, tmp2 = get_chain(case, resolution[1])
optimal_value = chain[p_index, argmin(tmp1)]
standard_deviation = std(chain[p_index,:])
minimum_value = quantile(chain[p_index, :], 0.05)
maximum_value = quantile(chain[p_index, :], 0.95)
range_value = (minimum_value, maximum_value)
if (optimal_value < maximum_value) && (optimal_value > minimum_value)
    #println("pass")
else
    println("-------")
    println(optimal_value)
    println(range_value)
    println("fail")
    println(case)
    println("-------")
end

end
###

scatter([1],[1], yerror = [(1,-1/2)])


###
case_range = 3:34 #34
mega_chain = randn(4, 32* 10001)
let tmp_chain = randn(4, 10001)
for resolution in resolutions[1:1]
    chains = []
    optimal_values = []
    standard_deviations = []
    stratifications = []
    range_values = []
    for case in cases[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if case == cases[case_range[1]]
            println("hi")
            tmp_chain = copy(chain)
        else
            println("augment")
            tmp_chain = combine(tmp_chain, chain)
        end
    end
end
@. mega_chain = tmp_chain
end
###
p = marginal_pdfs(mega_chain, left_bounds, right_bounds, parameter_dictionary, bins = 50)
p1 = plot(p...)

names = ["Surface Layer Fraction", "Nonlocal Amplitude", "Diffusivity Amplitude", "Unresolved Shear"]
left_bounds_j =  [0.0, 0.0, 0.0 , 4.0 ]
right_bounds_j = [0.1, 8.0, 2.0, 10.0]
p = joint_pdfs(mega_chain, left_bounds_j, right_bounds_j, parameter_dictionary, bins = 100)
p1 = plot(p[[5,3,6,2]]...)
median(mega_chain, dims = 2)
