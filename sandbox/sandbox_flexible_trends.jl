# use PyPlot backend
pyplot()

save_figures = true
###
# for stratification
bflux = L"[m/s^3]"
plot()
p_index = 5
range_list = [(0,0.02), (3,5),(1,2), (0,1), (0,1)]
case_range = 3:1:32   # 1:2 and 3:34. 32:34 isn't finished

y_range = range_list[p_index]
std_amplitude = 1
confidence_interval = 0.95
save_figures = true
rescale_p = true
Cᴿ = 0.3
if rescale_p
    range_list = [(0,0.02), (3,5),(1,2), (0,1), (0.0, 1 * Cᴿ) ]
end
y_range = range_list[p_index]
# loop over cases
for resolution in resolutions[1:1]
    chains = []
    optimal_values = []
    median_values = []
    standard_deviations = []
    stratifications = []
    range_values = []
    for case in cases[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        stratification = les.α * les.g * les.bottom_T
        push!(stratifications, stratification)
        filename = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
        mcmc_data = jldopen(filename, "r")
        chain = mcmc_data["𝑪"]
        if rescale_p
            @. chain[p_index, : ] *= Cᴿ
        end
        e1 = mcmc_data["ε"]
        e2 = mcmc_data["proposal_ε"]
        close(mcmc_data)
        push!(chains, chain)
        optimal_value = chain[p_index, argmin(e1)]
        println("case = $case")
        println(chain[:, argmin(e1)])
        median_value = median(chain[p_index,:])
        push!(optimal_values, optimal_value)
        push!(median_values, median_value)
        standard_deviation = std(chain[p_index,:])
        push!(standard_deviations, standard_deviation)

        minimum_value = median_value - quantile(chain[p_index, :], 1-confidence_interval)
        maximum_value = quantile(chain[p_index, :], confidence_interval) - median_value
        range_value = (minimum_value, maximum_value)
        push!(range_values, range_value)
    end
    N² = stratifications[:]
    C = optimal_values[:]
    C2 = median_values[:]
    σC = std_amplitude .* standard_deviations[:]
    res_string = @sprintf("%.2f ", 100/resolution[1])
    per_string = @sprintf("%.0f", confidence_interval* 100)
    p1 = scatter!(N², C2, xlabel = "Background Stratification, N² " * itime, ylabel = "Entrainment Coefficient", legend = :topleft, yerror = range_values, ylims = y_range, label = "Median values at " * res_string * "meter resolution")
    p1  = scatter!(N², C, label = "Optimal values at " * res_string * "meter resolution", shape = :star5, legend = :topleft, title = "Modes, Medians, and " * per_string * "% Probability Intervals" * " for h²N² scaling", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    if save_figures
        savefig(p1, pwd() * "/figures/new_scaling_trends.png")
    end
end

###
# for surface forcing
#flexible Surface
save_figures = true
plot()
p_index = 5
range_list = [(0,0.02), (3,5),(1,2), (0,1), (0,1)]
case_range = 1:10   # 1:2 and 3:34
#=
cr = []
push!(cr, collect(4:8:34))
push!(cr, collect(3:8:34))
push!(cr, collect(5:8:34))
case_range = [cr[1]; cr[2][2:end]; cr[3]]
=#

y_range = range_list[p_index]
std_amplitude = 1
confidence_interval = 0.95

# loop over cases
for resolution in resolutions[1:1]
    chains = []
    optimal_values = []
    median_values = []
    standard_deviations = []
    stratifications = []
    range_values = []
    for case in cases2[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        stratification =  les.top_T
        push!(stratifications, stratification)
        filename = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
        mcmc_data = jldopen(filename, "r")
        chain = mcmc_data["𝑪"]
        e1 = mcmc_data["ε"]
        e2 = mcmc_data["proposal_ε"]
        close(mcmc_data)
        push!(chains, chain)
        optimal_value = chain[p_index, argmin(e1)]
        println("case = $case")
        println(chain[:, argmin(e1)])
        median_value = median(chain[p_index,:])
        push!(optimal_values, optimal_value)
        push!(median_values, median_value)
        standard_deviation = std(chain[p_index,:])
        push!(standard_deviations, standard_deviation)

        minimum_value = median_value - quantile(chain[p_index, :], 1-confidence_interval)
        maximum_value = quantile(chain[p_index, :], confidence_interval) - median_value
        range_value = (minimum_value, maximum_value)
        push!(range_values, range_value)
    end
    N² = stratifications[:]
    C = optimal_values[:]
    C2 = median_values[:]
    σC = std_amplitude .* standard_deviations[:]
    res_string = @sprintf("%.2f ", 100/resolution[1])
    per_string = @sprintf("%.0f", confidence_interval* 100)
    p1 = scatter!(N², C2, xlabel = "Surface Buoyancy Flux " * bflux, ylabel = "Unresolved Shear", legend = :topleft, yerror = range_values, ylims = y_range, label = "Median values at " * res_string * "meter resolution")
    p1  = scatter!(N², C, label = "Optimal values at " * res_string * "meter resolution", shape = :star5, legend = :topleft, title = "Modes, Medians, and " * per_string * "% Confidence Intervals" * " for h²N² scaling", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    if save_figures
        savefig(p1, pwd() * "/figures/h2n2_surface.png")
    end
end
