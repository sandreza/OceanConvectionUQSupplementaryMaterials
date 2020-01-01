include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

# pdf figures
save_figures = false
case = cases[1]
#case = "compromise"
resolution = resolutions[1]

# get chain
chain, e1, e2 = get_chain(case, resolution[1])

# marginal pdfs
m,n = size(chain)
p = marginal_pdfs(chain[:,1:(n-1)], left_bounds, right_bounds, parameter_dictionary, bins = 100)
p1 = plot(p...)

if save_figures == true
    savefig(p1, pwd() * "/figures/figure_3.png")
end

# joint pdfs, reminders of bounds
names = ["Surface Layer Fraction", "Nonlocal Amplitude", "Diffusivity Amplitude", "Unresolved Shear"]
prob = 0.95
left_bounds_j = a_quantile(chain, 1-prob)
right_bounds_j = a_quantile(chain, prob)
p = joint_pdfs(chain, left_bounds_j, right_bounds_j, parameter_dictionary, bins = 150)
p1 = plot(p[[5,3,6,2]]...)

if save_figures == true
    savefig(p1, pwd() * "/figures/figure_3b.png")
end

###
plot(chain[4,:])
# extra
acceptance_rate = sum(e1 .== e2) / length(e1)
indmin = argmin(e1)
indmax = argmax(e1)
optimal_𝑪 = chain[:,indmin]

# means and standard deviations
μ, σ = calculate_partial_statistics(chain[:,1:10^6])
# index label
ind = 3
pμ, pσ = plot_partial_statistics(μ, σ; ind = ind);
plot(pμ)
plot(pσ)
# error in time and histogram
plot(e1)
histogram(e1[1:(10^6)], xlabel = "error", bins = 50, legend = false, normalize = true, ylabel = "pdf")
