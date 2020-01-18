include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

# compromise functions
# use PyPlot backend
pyplot()

save_figures = true

resolution = resolutions[1]
# define things for forward map
N = resolution[1]
Δt = resolution[2]
zᵖ = zeros(N)

ind_case_1 = 6
ind_case_2 = 34

# get LES 1
case = cases[ind_case_1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les1 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les1.t)
𝒢1 = CoreFunctionality.closure_free_convection(N, Δt, les1, subsample = subsample, grid = zᵖ)
ℒ1 = closure_default_loss_function(filename, N = N, Δt = Δt)
chain1, tmp1, tmp2 = get_chain(case, resolution[1])
indmin1 = argmin(tmp1)
𝑪1 = chain1[:,indmin1]
𝑪1 = median(chain1, dims = 2)[:]
println("median 1 is $𝑪1")
case = cases[ind_case_2]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les2 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les2.t)
𝒢2 = CoreFunctionality.closure_free_convection(N, Δt, les2, subsample = subsample, grid = zᵖ)
ℒ2 = closure_default_loss_function(filename, N = N, Δt = Δt)
chain2, tmp1, tmp2 = get_chain(case, resolution[1])
indmin2 = argmin(tmp1)
𝑪2 = chain2[:,indmin2]
𝑪2 = median(chain2, dims = 2)[:]
println("median 2 is $𝑪2")
# create other compromise distribution
chain4 = combine(chain1, chain2)
𝑪4 = median(chain4,dims=2)[:] #other choice of compromise


###
case_range1 = 1:(10^5-1)
case_range2 = copy(case_range1)
case_range4 = (5 * 10^4 +2) : (1 * 10^5 - 1 + 5 * 10^4)
index = 4
histogram(chain1[index,case_range1], normalize = true, alpha = 0.4, xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = 50, legend = true, ylabel = "pdf", label = "Case 1")
histogram!(chain2[index,case_range2], normalize = true, alpha = 0.4, label = "Case 2")
histogram!(chain4[index,case_range4], normalize = true, alpha = 0.4, label = "Compromise")
###
p = marginal_pdfs(chain1[:,case_range1], chain2[:,case_range2], left_bounds, right_bounds, parameter_dictionary)
plot(p...)
p1 = plot(p[4])

if save_figures == true
    savefig(p1, pwd() * "/figures/figure_9_distributions.pdf")
end
###
# parameters to loop over
# show everything in case 1 scenario
p_case1 = []
labels = ["Median 1", "Median 2",  "Compromise"]
parameter_list =[𝑪1, 𝑪2, 𝑪4]
maxT = maximum(les1.T[:,1])
minT = minimum(les1.T[:,1])
xlims = (minT,maxT)
plot()
for j in 1:3
    𝑪 = parameter_list[j]
    loss = ℒ1(𝑪)
    Tᵖ = 𝒢1(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les1.T[:,end], les1.z, label = "LES 1", legend = :bottomright, ylabel = "depth [m]", xlabel = "Temperature " * celsius, ylims = (-100, 0), xlims = xlims)
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 1 = " * loss_string * " " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case1,p1)
end
p1 = plot(p_case1...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_9_case_1.pdf")
end
###
# show everything in case 2 scenario
p_case2 = []
labels = ["Median 1", "Median 2", "Compromise"]
parameter_list =[𝑪1, 𝑪2, 𝑪4]
maxT = maximum(les2.T[:,1])
minT = minimum(les2.T[:,1])
xlims = (minT,maxT)
plot()
for j in 1:3
    𝑪 = parameter_list[j]
    loss = ℒ2(𝑪)
    Tᵖ = 𝒢2(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les2.T[:,end], les2.z, label = "LES 2", legend = :bottomright, ylabel = "depth [m]", xlabel = "Temperature " * celsius, ylims = (-100,0), xlims = xlims)
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 2 = " * loss_string * " " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case2,p1)
end
p2 = plot(p_case2...)
if save_figures == true
    savefig(p2, pwd() * "/figures/figure_9_case_2.pdf")
end
###
# Show median values for optimal in both cases
p1 = plot(p_case1[1], p_case2[2])
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_9_alternative.pdf")
end

# Show median values in opposite cases
p2 = plot(p_case1[2], p_case2[1])
if save_figures == true
    savefig(p2, pwd() * "/figures/figure_9.pdf")
end

# Show them together
p3 = plot(p_case1[1], p_case2[2], p_case1[2], p_case2[1])
if save_figures == true
    savefig(p3, pwd() * "/figures/figure_9_alternative3.pdf")
end
tmp = p3
