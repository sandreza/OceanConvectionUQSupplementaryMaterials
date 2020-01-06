include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

# compromise functions

save_figures = false
use_median = true
resolution = resolutions[1]
# define things for forward map
N = resolution[1]
Δt = resolution[2]
zᵖ = zeros(N)

ind_case_1 = 1
ind_case_2 = 2

# get LES 1
case = cases[ind_case_1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les1 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les1.t)
𝒢1 = CoreFunctionality.closure_free_convection_flexible(N, Δt, les1, subsample = subsample, grid = zᵖ, power = 1.0)
ℒ1 = CoreFunctionality.closure_T_nll(𝒢1, les1; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )
filename = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
indmin1 = argmin(e1)
chain1 = mcmc_data["𝑪"]
𝑪1 = chain1[:,indmin1]
if use_median
    𝑪1 = median(chain1, dims = 2)[:]
    println("median 1 is $𝑪1")
end
### now fo rcase 2
case = cases[ind_case_2]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les2 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les2.t)
𝒢2 = CoreFunctionality.closure_free_convection_flexible(N, Δt, les2, subsample = subsample, grid = zᵖ, power = 1.0)
ℒ2 = CoreFunctionality.closure_T_nll(𝒢2, les2; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )
filename = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
indmin2 = argmin(e1)
chain2 = mcmc_data["𝑪"]
𝑪2 = chain2[:,indmin2]
if use_median
    𝑪2 = median(chain2, dims = 2)[:]
    println("median 2 is $𝑪2")
end
# create other compromise distribution
chain4 = combine(chain1, chain2)
𝑪4 = median(chain4,dims=2)[:] #other choice of compromise

NN1 = 𝑪1[6]
NN2 = 𝑪2[6]

###
#=
case_range1 = 1:(10^3-1)
case_range2 = copy(case_range1)
case_range4 = (5 * 10^2 +2) : (1 * 10^3 - 1 + 5 * 10^2)
index = 5
histogram(chain1[index,case_range1], normalize = true, alpha = 0.4, xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = 50, legend = true, ylabel = "pdf", label = "Case 1")
histogram!(chain2[index,case_range2], normalize = true, alpha = 0.4, label = "Case 2")
histogram!(chain4[index,case_range4], normalize = true, alpha = 0.4, label = "Compromise")
###
p = marginal_pdfs(chain1[:,case_range1], chain2[:,case_range2], left_bounds, right_bounds, parameter_dictionary)
plot(p...)
p1 = plot(p[4])

if save_figures == true
    savefig(p1, pwd() * "/figures/figure_9_distributions.png")
end
=#
###
# parameters to loop over
# show everything in case 1 scenario
p_case1 = []
labels = ["Mode 1", "Mode 2",  "Compromise"]
if use_median
    labels = ["Median 1", "Median 2",  "Compromise"]
end
#𝑪1[1] = 1e-4
#𝑪2[2] = 3.5
#𝑪2[5] = 0.375
#𝑪2[3] = 10
parameter_list =[𝑪1, 𝑪2, 𝑪4]
# the stratification changes
𝑪1[6] = NN1
𝑪2[6] = NN1
𝑪4[6] = NN1
maxT = maximum(les1.T[:,end] .+ 0.1)
minT = minimum(les1.T[:,end])
plot()
for j in 1:3
    𝑪 = parameter_list[j]
    loss = maximum( ℒ1(𝑪) )
    Tᵖ = 𝒢1(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les1.T[:,end], les1.z, label = "LES 1", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature [C]", ylims = (-100, 0), xlims = (minT, maxT))
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 1 = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case1,p1)
end
p1 = plot(p_case1...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_9_case_1.png")
end
###
# show everything in case 2 scenario
p_case2 = []
labels = ["Mode 1", "Mode 2",  "Compromise"]
if use_median
    labels = ["Median 1", "Median 2",  "Compromise"]
end
parameter_list =[𝑪1, 𝑪2, 𝑪4]
𝑪1[6] = NN2
𝑪2[6] = NN2
𝑪4[6] = NN2
plot()
maxT = maximum(les2.T[:,end] .+ 0.1)
minT = minimum(les2.T[:,end])
xlims = (minT,maxT)
for j in 1:3
    𝑪 = parameter_list[j]
    loss = maximum(ℒ2(𝑪))
    Tᵖ = 𝒢2(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les2.T[:,end], les2.z, label = "LES 2", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature [C]", ylims = (-100,0), xlims = (minT,maxT))
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 2 = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case2,p1)
end
p2 = plot(p_case2...)
if save_figures == true
    savefig(p2, pwd() * "/figures/figure_11_case_2.png")
end
###
# Show median values for optimal in both cases
p1 = plot(p_case1[1], p_case2[2])
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_11_alternative.png")
end

# Show median values in opposite cases
p2 = plot(p_case1[2], p_case2[1])
if save_figures == true
    savefig(p2, pwd() * "/figures/figure_11.png")
end

# Show them together
p4 = plot(p_case1[1], p_case2[2], p_case1[2], p_case2[1])
if save_figures == true
    savefig(p4, pwd() * "/figures/figure_11_alternative3.png")
end

# Show compromise values in both cases
p5 = plot(p_case1[3], p_case2[3])
if save_figures == true
    savefig(p5, pwd() * "/figures/figure_11_alternative4.png")
end
