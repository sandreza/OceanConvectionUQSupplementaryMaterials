include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

# use PyPlot backend
pyplot()



save_figures = true
rescale_4 = true
Cᴿ = 0.3
#trends plots

case_range = 3:34
y_range = (0,16)
if rescale_4
    y_range = (0, 16* Cᴿ)
end
std_amplitude = 1
chains = []
optimal_values = []
standard_deviations = []
stratifications = []

plot()
p_index = 4

for resolution in resolutions[1:3]
    chains = []
    optimal_values = []
    standard_deviations = []
    stratifications = []

    for case in cases[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        stratification = les.α * les.g * les.bottom_T
        push!(stratifications, stratification)
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if rescale_4
            chain[4,:] *= Cᴿ
        end
        push!(chains, chain)
        optimal_value = chain[p_index, argmin(tmp1)]
        push!(optimal_values, optimal_value)
        standard_deviation = std(chain[p_index,:])
        push!(standard_deviations, standard_deviation)
    end
    N² = stratifications[:]
    C = optimal_values[:]
    σC = std_amplitude .* standard_deviations[:]
    res_string = @sprintf("%.2f ", 100/resolution[1])
    p1 = scatter!(N², C, xlabel = "Background Stratification, N^2 " * itime, ylabel = names[p_index], legend = :topleft, ribbon = σC, fillalpha = 0.2, ylims = y_range, label = "Resolution = " * res_string * "meters")
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_7.png")
    end
end

###
plot()
p_index = 4

case_range = 3:34
y_range = (0, 16.0)
if rescale_4
    y_range = (0, 16* Cᴿ)
end
std_amplitude = 1
for resolution in resolutions[1:3]
    chains = []
    optimal_values = []
    standard_deviations = []
    stratifications = []
    range_values = []
    for case in cases[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        stratification = les.α * les.g * les.bottom_T
        push!(stratifications, stratification)
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if rescale_4
            chain[4,:] *= Cᴿ
        end
        push!(chains, chain)
        optimal_value = chain[p_index, argmin(tmp1)]
        push!(optimal_values, optimal_value)
        standard_deviation = std(chain[p_index,:])
        push!(standard_deviations, standard_deviation)
        maximum_value = maximum(chain[p_index, :]) - optimal_value
        minimum_value = optimal_value - minimum(chain[p_index, :])
        range_value = (minimum_value, maximum_value)
        push!(range_values, range_value)
    end
    N² = stratifications[:]
    println(N²)
    ΔN² = (maximum(N²)-minimum(N²))/10
    minN² = minimum(N²)
    maxN² = maximum(N²)
    ticks = maxN²:ΔN²:minN²
    C = optimal_values[:]
    σC = std_amplitude .* standard_deviations[:]
    res_string = @sprintf("%.2f ", 100/resolution[1])
    p1 = scatter!(N², C, xlabel = "Background Stratification, N^2 " * itime, ylabel = names[p_index], legend = :topleft, yerror = σC, ylims = y_range, label = "Resolution = " * res_string * "meters")
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_7_alternative.png")
    end
end

###

#
#
plot()
p_index = 4

case_range = 3:34
y_range = (0,16)
if rescale_4
    y_range = (0, 16* Cᴿ)
end
std_amplitude = 1
confidence_interval = 0.95
for resolution in resolutions[1:1]
    chains = []
    optimal_values = []
    standard_deviations = []
    stratifications = []
    range_values = []
    for case in cases[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        stratification = les.α * les.g * les.bottom_T
        push!(stratifications, stratification)
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if rescale_4
            chain[4,:] *= Cᴿ
        end
        push!(chains, chain)
        #optimal_value = chain[p_index, argmin(tmp1)]
        optimal_value = median(chain[p_index,:])
        push!(optimal_values, optimal_value)
        standard_deviation = std(chain[p_index,:])
        push!(standard_deviations, standard_deviation)
        minimum_value = optimal_value - quantile(chain[p_index, :], 1-confidence_interval)
        maximum_value = quantile(chain[p_index, :], confidence_interval) - optimal_value
        range_value = (minimum_value, maximum_value)
        push!(range_values, range_value)
    end
    N² = stratifications[:]
    C = optimal_values[:]
    σC = std_amplitude .* standard_deviations[:]
    res_string = @sprintf("%.2f ", 100/resolution[1])
    p1 = scatter!(N², C, xlabel = "Background Stratification, N^2 " * itime, ylabel = names[p_index], legend = :topleft, yerror = range_values, ylims = y_range, label = "Resolution = " * res_string * "meters")
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_7_alternative2.png")
    end
end
###
plot()
p_index = 4
range_list = [(0,1), (0,8),(0,12), (0,16)]
if rescale_4
    range_list = [(0,1), (0,8),(0,12), (0, 16* Cᴿ)]
end
case_range = 3:34   # 1:2 and 3:34
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
    for case in cases[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        stratification = les.α * les.g * les.bottom_T
        push!(stratifications, stratification)
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if rescale_4
            chain[4,:] *= Cᴿ
        end
        push!(chains, chain)
        optimal_value = chain[p_index, argmin(tmp1)]
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
    per_string = @sprintf("%.0f", (2*confidence_interval-1)* 100)
    p1 = scatter!(N², C2, xlabel = "Background Stratification, N² " * itime, ylabel = names[p_index], legend = :topleft, yerror = range_values, ylims = y_range, label = "Median values at " * res_string * "meter resolution")
    p1  = scatter!(N², C, label = "Optimal values at " * res_string * "meter resolution", shape = :star5, legend = :topleft, title = "Modes, Medians, and " * per_string * "% Probability Intervals", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_7_alternative3.png")
    end
end


###

# include all together


p_index = 4
range_list = [(0,1), (0,8),(0,6), (0,16)]
if rescale_4
    range_list = [(0,1), (0,8),(0,6), (0,16*Cᴿ)]
end
case_range = 3:34   # 1:2 and 3:34

std_amplitude = 1
confidence_interval = 0.95
p = []
# loop over cases
for p_index in [1,2,3,4]
    plot()
    y_range = range_list[p_index]
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
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if rescale_4
            chain[4,:] *= Cᴿ
        end
        push!(chains, chain)
        optimal_value = chain[p_index, argmin(tmp1)]
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
    p1 = scatter!(N², C2, xlabel = "Background Stratification, N² " * itime, ylabel = names[p_index], legend = false, yerror = range_values, ylims = y_range, label = "Median values at " * res_string * "meter resolution")
    p1  = scatter!(N², C, label = "Optimal values at " * res_string * "meter resolution", shape = :star5, legend = false, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_7_alternative_" * names[p_index] * ".png")
    end
    push!(p,p1)
end
p1 = plot(p...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_7_alternative5.png")
end
