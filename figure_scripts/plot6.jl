include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

# compromise functions

save_figures = true

resolution = resolutions[1]
# define things for forward map
N = resolution[1]
Δt = resolution[2]
zᵖ = zeros(N)

# get LES 1
case = cases[1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les1 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les1.t)
𝒢1 = CoreFunctionality.closure_free_convection(N, Δt, les1, subsample = subsample, grid = zᵖ)
ℒ1 = closure_default_loss_function(filename, N = N, Δt = Δt)
chain1, tmp1, tmp2 = get_chain(case, resolution[1])
indmin1 = argmin(tmp1)
𝑪1 = chain1[:,indmin1]
# get LES 2
case = cases[2]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les2 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les2.t)
𝒢2 = CoreFunctionality.closure_free_convection(N, Δt, les2, subsample = subsample, grid = zᵖ)
ℒ2 = closure_default_loss_function(filename, N = N, Δt = Δt)
chain2, tmp1, tmp2 = get_chain(case, resolution[1])
indmin2 = argmin(tmp1)
𝑪2 = chain2[:,indmin2]

# get compromise data
case = "compromise"
chain3, e1, e2 = get_chain(case, resolution[1])
indmin3 = argmin(e1)
𝑪3 = chain3[:,indmin3]

# create other compromise distribution
chain4 = combine(chain1, chain2)
𝑪4 = median(chain4,dims=2)[:] #other choice of compromise

# construct compromise loss function
a = ℒ1(𝑪1)
b = ℒ1(𝑪2)
c = ℒ2(𝑪2)
d = ℒ2(𝑪1)
# now define combined loss function
scale = (c+d) / (a+b)
ℒ_compromise(𝑪) = 0.5 *( ℒ1(𝑪) * scale + ℒ2(𝑪) )

###
case_range1 = 1:(10^6-1)
case_range2 = copy(case_range1)
case_range3 = 1:(10^6-1)
case_range4 = (5 * 10^5 +2) : (1 * 10^6 - 1 + 5 * 10^5)
index = 4
histogram(chain1[index,case_range1], normalize = true, alpha = 0.4, xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = 50, legend = true, ylabel = "pdf", label = "Case 1")
histogram!(chain2[index,case_range2], normalize = true, alpha = 0.4, label = "Case 2")
histogram!(chain3[index,case_range3], normalize = true, alpha = 0.4, label = "Compromise 1")
histogram!(chain4[index,case_range4], normalize = true, alpha = 0.4, label = "Compromise 2")
###
p = marginal_pdfs(chain1[:,case_range1], chain2[:,case_range2], left_bounds, right_bounds, parameter_dictionary)
plot(p...)
p1 = plot(p[4])

if save_figures == true
    savefig(p1, pwd() * "/figures/figure_6_distributions.png")
end
###
ℒ_compromise(𝑪3)

###
# parameters to loop over
# show everything in case 1 scenario
p_case1 = []
labels = ["Optimal 1", "Optimal 2", "Compromise", "Compromise 2"]
parameter_list =[𝑪1, 𝑪2, 𝑪3, 𝑪4]
plot()
for j in 1:4
    𝑪 = parameter_list[j]
    loss = ℒ1(𝑪)
    Tᵖ = 𝒢1(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les1.T[:,end], les1.z, label = "LES 1", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature [C]", ylims = (-90, 0), xlims = (19.1, 19.45))
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 1 = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case1,p1)
end
p1 = plot(p_case1...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_6_case_1.png")
end
###
# show everything in case 2 scenario
p_case2 = []
labels = ["Optimal 1", "Optimal 2", "Compromise", "Compromise 2"]
parameter_list =[𝑪1, 𝑪2, 𝑪3, 𝑪4]
plot()
for j in 1:4
    𝑪 = parameter_list[j]
    loss = ℒ2(𝑪)
    Tᵖ = 𝒢2(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les2.T[:,end], les2.z, label = "LES 2", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature [C]", ylims = (-60,0), xlims = (17.5, 18.5))
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 2 = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case2,p1)
end
p2 = plot(p_case2...)
if save_figures == true
    savefig(p2, pwd() * "/figures/figure_6_case_2.png")
end
###
# show both cases at once with compromise optimal
p1 = plot(p_case1[3], p_case2[3])
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_6.png")
end

p1 = plot(p_case1[1], p_case2[2], p_case1[3], p_case2[3])
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_6_alternative.png")
end


p1 = plot(p_case1[1], p_case2[2], p_case1[2], p_case2[1])
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_6_alternative_2.png")
end
