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
Cᴿ = 0.3
chain[4,:] *= Cᴿ #mutliply by the critical richardson number, for visualizing pdfs
right_bounds[4] *= Cᴿ

# marginal pdfs
m,n = size(chain)
p = marginal_pdfs(chain[:,1:(n-1)], left_bounds, right_bounds, parameter_dictionary, bins = 50)
p1 = plot(p...)

if save_figures == true
    savefig(p1, pwd() * "/figures/figure_3.png")
end

# joint pdfs, reminders of bounds
prob = 0.95
left_bounds_j = a_quantile(chain, 1-prob)
right_bounds_j = a_quantile(chain, prob)
p = joint_pdfs(chain, left_bounds_j, right_bounds_j, parameter_dictionary, bins = 50)
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
μ, σ = calculate_partial_statistics(chain[:,1:end])
# index label
ind = 4
pμ, pσ = plot_partial_statistics(μ, σ; ind = ind);
plot(pμ)
plot(pσ)
# error in time and histogram
plot(e1)
histogram(e1[1:(end-1)], xlabel = "error", bins = 50, legend = false, normalize = true, ylabel = "pdf")
