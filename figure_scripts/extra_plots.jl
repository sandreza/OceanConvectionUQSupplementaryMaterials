include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")
pyplot()
parameter_list = [default_𝑪, optimal_𝑪, mean_𝑪, median_𝑪]
plot()
p = []
xup = 19.45
xdown = 19.2
yup = 0.0
ydown = -75
Δy = abs(yup-ydown)
Δx = abs(xup -xdown)

p = []

# get the les data
case = cases[1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)
subsample_parameter = 1
start = 1
subsample = start:subsample_parameter:length(les.t)
##
N = Nlist[1]
Δt = les.t[2] - les.t[1]
zᵖ = zeros(N)
# define the forward map
𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample, grid = zᵖ)
# define the loss function
ℒᵗ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )


p1 = plot(les.T[:,end], les.z, label = "LES", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius, xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)

j = 1
𝑪 = parameter_list[j]
loss_default = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
daystring = @sprintf("%.1f", les.t[end] ./ 86400)
p1 = scatter!(Tᵖ[:,end], zᵖ, label = "Default", title = "t = "* daystring * " days", markersize = 6, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokealpha = 0.0,  markerstrokewidth = 0.0, markershape = :square, aspect_ratio = 2 * Δx/Δy, markercolor = :green)

j = 2
𝑪 = parameter_list[j]
loss_optimal = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
p1 = scatter!(Tᵖ[:,end], zᵖ, label = "Mode", title = "t = "* daystring * " days", markersize = 6, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokealpha = 0.0,  markerstrokewidth = 0.0, markershape = :circle, aspect_ratio = 2 * Δx/Δy, markercolor = :magenta)

inds = 30:length(les.t)
p2 = plot(les.t[inds] ./ 86400, sqrt.(loss_default[inds]), legend = :bottomright, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Default" , color = :green)
p2 =  plot!(les.t[inds] ./ 86400, sqrt.(loss_optimal[inds]), legend = :topleft, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Mode" , color = :magenta, title = "Error in Time")


p = plot(p1, p2)
savefig(p, pwd() * "/figures/profile_and_loss.pdf")

###

# include("plot2.jl")
plot()
anim = @animate for i in 1:50:length(les.t)
    xup = 20.0 #maximum(les.T[:,i]) + 0.05 # 20.0 - les.t[i] /( 8 *86400) * 0.55
    xdown = 19.2 #19.65 - les.t[i] /( 8 *86400) * 0.45
    yup = 0.0
    ydown = -75 # -30 - les.t[i] /( 8 * 86400) * 45
    Δy = 75.0 # abs(yup-ydown)
    Δx = 0.25 # abs(xup -xdown)
    Δy = abs(yup-ydown)
    Δx = abs(xup -xdown)
    daystring = @sprintf("%.0f", les.t[end] ./ 86400)
    p1 = plot(les.T[:,i], les.z, label = "LES", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius, xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, title = "t = "* daystring * " days")
    j = 1
    𝑪 = parameter_list[j]
    loss_default = ℒᵗ(𝑪)
    Tᵖ = 𝒢(𝑪)

    p2 = scatter(Tᵖ[:,i], zᵖ, label = labels[j], title = "t = "* daystring * " days", markersize = 6, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokealpha = 0.0,  markerstrokewidth = 0.0, markershape = :circle, aspect_ratio = 2 * Δx/Δy, markercolor = :blue, xlims = (xdown, xup), ylims = (ydown, yup), ylabel = "depth [m]", xlabel = "Temperature " * celsius, legend = :topleft)

    j = 2
    𝑪 = parameter_list[j]
    loss_default = ℒᵗ(𝑪)
    Tᵖ = 𝒢(𝑪)
    p3 = scatter(Tᵖ[:,i], zᵖ, label = labels[j], title = "t = "* daystring * " days", markersize = 6, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokealpha = 0.0,  markerstrokewidth = 0.0, markershape = :square, aspect_ratio = 2 * Δx/Δy, markercolor = :green, xlims = (xdown, xup), ylims = (ydown, yup), ylabel = "depth [m]", xlabel = "Temperature " * celsius, legend = :topleft)

    pp = plot(p2, p3, p1, layout = (1,3))
# display(pp)
end

gif(anim, pwd() * "/figures/ocean_sciences_dynamic.gif", fps = 15)
###
ti = argmin(les.t .< 2*86400)
les.t[ti]./86400
tfield = les.T[:,ti]
scale = 10^2
bfield = @. tfield * les.α * les.g * scale
xup = 19.725
xdown = 19.5
yup = 0.0
ydown = -50

Δy = abs(yup-ydown)
Δx = abs(xup -xdown)

plot(tfield, les.z, ylims = (ydown,yup), xlims = (xdown, xup), linewidth = 2 , grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box,  ylabel = "Depth [m]", xlabel = "Temperature " * celsius, legend = false, title = " ", top_margin = 7.0mm, aspect_ratio = 2 * Δx/Δy)


btop = les.α * les.g * xup * scale
bbottom = les.α * les.g * xdown * scale

Δx *= les.α * les.g * scale
#p1 = plot!(twinx(), bfield, les.z, ylims = (bbottom, btop), xlims = (bbottom, btop), xlabel = "Buoyancy " * acceleration, legend = false, xmirror = true, yaxis = false, linewidth = 2 , grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = 2 * Δx/Δy)

p1 = plot!(twinx(), bfield, les.z, ylims = (ydown, yup), xlims = (bbottom, btop), xlabel = "Buoyancy x 10² " * acceleration, legend = false, xmirror = true, yaxis = false, linewidth = 2 , grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, aspect_ratio = 2 * Δx/Δy)

plot(p1)
savefig(p1, pwd() * "/figures/t_prof.pdf")

savefig(p1, pwd() * "/figures/buoyancy_and_temperature.pdf")
#=
tindex = argmin(abs.(les.t .- 86400 *2))
plot(les.T[:,tindex], les.z, label = "LES", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius, xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
=#

###
plot()
save_figures = true
rescale_4 = true
p_index = 4
Cᴿ = 0.3
range_list = [(0,1), (0,8),(0,12), (0,16)]
default_val = default_𝑪[p_index]
if rescale_4
    range_list = [(0,1), (0,8),(0,12), (0, 16* Cᴿ)]
    default_val = default_𝑪[p_index] * Cᴿ
end
case_range = 3:34   # 1:2 and 3:34
y_range = range_list[p_index]
std_amplitude = 1
confidence_interval = 0.95
names = [L"C^S", L"C^N", L"C^D", L"C^H", L"C^H"]
# loop over cases

resolution = resolutions[1]
chains = []
optimal_values = []
median_values = []
standard_deviations = []
stratifications = []
range_values = []
plot()
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

p1 = plot!(N², 0 * N² .+ default_val , color = :black, linewidth = 2, label = "Default", linestyle = :dash)

p1 = scatter!(N², C2, xlabel = "Background Stratification, N² " * itime, ylabel = names[p_index], yerror = range_values, ylims = y_range, color = :blue, markersize = 6, label = "Median ")

p1  = scatter!(N², C, label = "Mode ", shape = :star5, legend = :topleft, title = "Modes, Medians, and " * per_string * "% Probability Intervals", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokewidth = 0.0, markersize = 7,  markercolor = :red, markerstrokealpha = 0.0)
if rescale_4 && (p_index==4)
    default_val = default_𝑪[p_index] * Cᴿ
end
display(p1)
if save_figures == true
    savefig(p1, pwd() * "/figures/kpp_default_trends.pdf")
end

###



plot()
rescale_4 = true
p_index = 5
Cᴿ = 0.3
range_list = [(0,0.02), (3,5),(1,2), (0,1), (0,1)]
if rescale_4
    range_list = [(0,0.02), (3,5),(1,2), (0,1), (0,1*Cᴿ)]
end
names = ["Surface Layer Fraction", "Nonlocal Amplitude", "Diffusivity Amplitude", "Mixing Parameter", "Mixing Parameter", "Exponent"]

names = [L"C^S", L"C^N", L"C^D", L"C^\star", L"C^\star"]

case_range = 3:34   # 1:2 and 3:34
y_range = range_list[p_index]
std_amplitude = 1
confidence_interval = 0.95

# loop over cases
resolution = resolutions[1]
chains = []
optimal_values = []
median_values = []
standard_deviations = []
stratifications = []
range_values = []
plot()
for case in cases[case_range]
    filename = pwd() * "/LES/" * case * "_profiles.jld2"
    les = CoreFunctionality.OceananigansData(filename)
    stratification = les.α * les.g * les.bottom_T
    push!(stratifications, stratification)
    filename = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
    mcmc_data = jldopen(filename, "r")
    chain = mcmc_data["𝑪"]
    tmp1 = mcmc_data["ε"]
    tmp2 = mcmc_data["proposal_ε"]
    close(mcmc_data)
    if rescale_4
        chain[5,:] *= Cᴿ
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

p1 = plot!(N², 0 * N² .+ 1/6 , color = :black, linewidth = 2, label = "Theoretical ", linestyle = :dash)

p1 = scatter!(N², C2, xlabel = "Background Stratification, N² " * itime, ylabel = names[p_index], yerror = range_values, ylims = y_range, color = :blue, markersize = 6, label = "Median ")

p1  = scatter!(N², C, label = "Mode ", shape = :star5, legend = :topleft, title = "Modes, Medians, and " * per_string * "% Probability Intervals", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokewidth = 0.0, markersize = 7,  markercolor = :red, markerstrokealpha = 0.0)

display(p1)
if save_figures == true
    savefig(p1, pwd() * "/figures/kpp_new_trends.pdf")
end
