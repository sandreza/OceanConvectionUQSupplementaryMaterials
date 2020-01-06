include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")

save_figures = true
rescale_4 = true
Cᴿ = 0.3

case_range = 1:10
y_range = (1,16)
if rescale_4
    y_range = (0, 16* Cᴿ)
end
std_amplitude = 1
chains = []
optimal_values = []
standard_deviations = []
surface_forcings = []

plot()
p_index = 4

for resolution in resolutions[1:3]
    chains = []
    optimal_values = []
    standard_deviations = []
    surface_forcings = []

    for case in cases2[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        surface_forcing = les.α * les.g * les.top_T
        push!(surface_forcings, surface_forcing)
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
    Φ = surface_forcings[:]
    C = optimal_values[:]
    σC = std_amplitude .* standard_deviations[:]
    res_string = @sprintf("%.2f ", 100/resolution[1])
    p1 = scatter!(Φ, C, xlabel = "Bouyancy Forcing, s^(-2)", ylabel = names[p_index], legend = :topright, ribbon = σC, fillalpha = 0.2, ylims = y_range, label = "Resolution = " * res_string * "meters")
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_8.png")
    end
end
###


case_range = 1:10
y_range = (1,16)
if rescale_4
    y_range = (0, 16* Cᴿ)
end
std_amplitude = 1
chains = []
optimal_values = []
standard_deviations = []
surface_forcings = []

plot()
p_index = 4
confidence_interval = 0.95
for resolution in resolutions[1:1]
    chains = []
    optimal_values = []
    standard_deviations = []
    surface_forcings = []
    range_values = []
    median_values = []
    for case in cases2[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        surface_forcing = les.α * les.g * les.top_T
        push!(surface_forcings, surface_forcing)
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if rescale_4
            chain[4,:] *= Cᴿ
        end
        push!(chains, chain)
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
    Φ = surface_forcings[:]
    C = optimal_values[:]
    C2 = median_values[:]
    res_string = @sprintf("%.2f ", 100/resolution[1])
    per_string = @sprintf("%.0f", confidence_interval* 100)
    σC = std_amplitude .* standard_deviations[:]

    p1 = scatter!(Φ, C2, xlabel = "Surface Bouyancy Forcing, [1/s²]", ylabel = names[p_index], legend = :topleft, yerror = range_values, ylims = y_range, label = "Medians at " * res_string * "meter resolution")

    p1 = scatter!(Φ, C, xlabel = "Surface Bouyancy Forcing, [1/s²]", ylabel = names[p_index], label = "Modes at " * res_string * "meter resolution", shape = :star5, legend = :topright, title = "Modes, Medians, and " * per_string * "% Confidence Intervals", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_8_alternative.png")
    end
end


###
plot()
p_index = 4
range_list = [(0,1), (0,8),(0,12), (0,16)]
if rescale_4
    range_list = [(0,1), (0,8),(0,12), (0,16 * Cᴿ)]
end
case_range = 1:10
y_range = range_list[p_index]
std_amplitude = 1
probability = 0.99
for resolution in resolutions[1:1]
    chains = []
    optimal_values = []
    median_values = []
    standard_deviations = []
    surface_forcings = []
    range_values = []
    for case in cases2[case_range]
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        les = CoreFunctionality.OceananigansData(filename)
        surface_forcing = les.α * les.g * les.top_T
        push!(surface_forcings, surface_forcing)
        chain, tmp1, tmp2 = get_chain(case, resolution[1])
        if rescale_4
            chain[4,:] *= Cᴿ
        end
        push!(chains, chain)
        optimal_value = chain[p_index, argmin(tmp1)]
        median_value = median(chain[p_index,:])
        push!(optimal_values, optimal_value)
        push!(median_values,median_value)
        standard_deviation = std(chain[p_index,:])
        push!(standard_deviations, standard_deviation)
        minimum_value = median_value - quantile(chain[p_index, :], 1-probability )
        maximum_value =  quantile(chain[p_index, :], probability) -median_value
        range_value = (minimum_value, maximum_value)
        push!(range_values, range_value)
    end
    Φ = surface_forcings[:]
    C = median_values[:]
    C2 = optimal_values[:]
    σC = std_amplitude .* standard_deviations[:]
    res_string = @sprintf("%.2f ", 100.0/resolution[1])
    p1 = scatter!(Φ, C, xlabel = "Surface Bouyancy Forcing, s^(-2)", ylabel = names[p_index], legend = :topright, yerror = range_values, fillalpha = 0.2, ylims = y_range, label = "Medians at Resolution = " * res_string * "meters", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    p1  = scatter!(Φ, C2, label = "Modes at Resolution = " * res_string * "meters", shape = :star5)
    display(p1)
    if save_figures == true
        savefig(p1, pwd() * "/figures/figure_8_alternative2.png")
    end
end
