include(pwd() * "/src/LocalOceanUQSupplementaryMaterials.jl")
include(pwd() * "/scripts/utils.jl")
include(pwd() * "/figure_scripts/utils.jl")
pyplot()
##
case = cases[1]
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
# define the loss function
ℒ = CoreFunctionality.closure_T_nll(𝒢, les)
# define time dependent loss function
ℒᵗ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )

# get MCMC data
resolution_label = "_res_" * string(N)
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
acceptance_rate = sum(e1 .== e2) / length(e1)
println("the acceptance rate was")
println(acceptance_rate)
indmin = argmin(e1)
close(mcmc_data)

seconds_in_a_day = 86400
# parameters to loop over
labels = ["Default", "Mode", "Mean", "Median"]
default_𝑪 = [0.1, 6.33, 1.36, 3.19]
#default_𝑪 = [0.0760809666611145; 4.342473912404762; 2.1630355831002954; 5.57111619953263] # from megachain
optimal_𝑪 = chain[:, indmin]
mean_𝑪 = mean(chain,dims=2)
median_𝑪 = median(chain,dims=2)

##
parameter_list = [default_𝑪, optimal_𝑪, mean_𝑪, median_𝑪]

# get the les data
case = cases[1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)
subsample_parameter = 1
start = 1
subsample = start:subsample_parameter:length(les.t)
##
N = N
Δt = les.t[2] - les.t[1]
zᵖ = zeros(N)
# define the forward map
𝒢 = CoreFunctionality.closure_free_convection(N, Δt, les, subsample = subsample, grid = zᵖ)
# define the loss function
ℒᵗ = CoreFunctionality.closure_T_nll(𝒢, les; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )
##
anim = @animate for i in 30:10:length(les.t)
xup = maximum(les.T[:,i]) + 0.04
xdown = 19.2
yup = 0.0
ydown = -75
Δy = abs(yup-ydown)
Δx = abs(xup -xdown)

p1 = plot(les.T[:,i], les.z, label = "LES", legend = :topleft, 
ylabel = "depth [m]", xlabel = "Temperature " * celsius, 
xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, 
aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, 
gridalpha = 0.25, framestyle = :box)
j = 1
𝑪 = parameter_list[j]
loss_default = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
daystring = @sprintf("%.1f", les.t[i] ./ 86400)
p1 = scatter!(Tᵖ[:,i], zᵖ, label = "Reference", 
title = "t = "* daystring * " days", markersize = 6, 
grid = true, gridstyle = :dash, gridalpha = 0.25, 
framestyle = :box, markerstrokealpha = 0.0,  
markerstrokewidth = 0.0, markershape = :square,
 aspect_ratio = 2 * Δx/Δy, markercolor = :green)
j = 2
𝑪 = parameter_list[j]
loss_optimal = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
p1 = scatter!(Tᵖ[:,i], zᵖ, label = "Mode", title = "t = "* daystring * " days", markersize = 6, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokealpha = 0.0,  markerstrokewidth = 0.0, markershape = :circle, aspect_ratio = 2 * Δx/Δy, markercolor = :magenta)
inds = 30:i
p2 = plot(les.t[inds] ./ 86400, sqrt.(loss_default[inds]), legend = :bottomright, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Reference" , color = :green)
p2 =  plot!(les.t[inds] ./ 86400, sqrt.(loss_optimal[inds]), legend = :topleft, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Mode" , color = :magenta, title = "Error in Time")

p = plot(p1, p2)
end
gif(anim, pwd() * "/error_animation.gif", fps = 15)

##
anim = @animate for i in 30:10:length(les.t)
xup = maximum(les.T[:,i]) + 0.04
xup = 20.01
xdown = 19.2
yup = 0.0
ydown = -75
Δy = abs(yup-ydown)
Δx = abs(xup -xdown)

p1 = plot(les.T[:,i], les.z, label = "LES", legend = :bottomright, 
ylabel = "depth [m]", xlabel = "Temperature " * celsius, 
xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, 
aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, 
gridalpha = 0.25, framestyle = :box)
j = 1
𝑪 = parameter_list[j]
loss_default = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
daystring = @sprintf("%.1f", les.t[i] ./ 86400)
p1 = scatter!(Tᵖ[:,i], zᵖ, label = "Reference", 
title = "t = "* daystring * " days", markersize = 6, 
grid = true, gridstyle = :dash, gridalpha = 0.25, 
framestyle = :box, markerstrokealpha = 0.0,  
markerstrokewidth = 0.0, markershape = :square,
 aspect_ratio = 2 * Δx/Δy, markercolor = :green)
j = 2
𝑪 = parameter_list[j]
loss_optimal = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
p1 = scatter!(Tᵖ[:,i], zᵖ, label = "Mode", title = "t = "* daystring * " days", markersize = 6, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokealpha = 0.0,  markerstrokewidth = 0.0, markershape = :circle, aspect_ratio = 2 * Δx/Δy, markercolor = :magenta)
inds = 30:i
p2 = plot(les.t[inds] ./ 86400, sqrt.(loss_default[inds]), legend = :bottomright, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Reference" , color = :green, ylims = (-0.0001, 0.017), xlims = (0.0, 8.0))
p2 =  plot!(les.t[inds] ./ 86400, sqrt.(loss_optimal[inds]), legend = :topleft, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Mode" , color = :magenta, title = "Error in Time")

p = plot(p1, p2)
end
gif(anim, pwd() * "/error_animation_fixed.gif", fps = 15)
mp4(anim, pwd() * "/error_animation_fixed.mp4", fps = 15)

##
# without the optimum
anim = @animate for i in 30:10:length(les.t)
xup = maximum(les.T[:,i]) + 0.04
xup = 20.01
xdown = 19.2
yup = 0.0
ydown = -75
Δy = abs(yup-ydown)
Δx = abs(xup -xdown)

p1 = plot(les.T[:,i], les.z, label = "LES", legend = :bottomright, 
ylabel = "depth [m]", xlabel = "Temperature " * celsius, 
xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, 
aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, 
gridalpha = 0.25, framestyle = :box)
j = 1
𝑪 = parameter_list[j]
loss_default = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
daystring = @sprintf("%.1f", les.t[i] ./ 86400)
p1 = scatter!(Tᵖ[:,i], zᵖ, label = "Reference", 
title = "t = "* daystring * " days", markersize = 6, 
grid = true, gridstyle = :dash, gridalpha = 0.25, 
framestyle = :box, markerstrokealpha = 0.0,  
markerstrokewidth = 0.0, markershape = :square,
 aspect_ratio = 2 * Δx/Δy, markercolor = :green)
j = 2
𝑪 = parameter_list[j]
loss_optimal = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
# p1 = scatter!(Tᵖ[:,i], zᵖ, label = "Mode", title = "t = "* daystring * " days", markersize = 6, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, markerstrokealpha = 0.0,  markerstrokewidth = 0.0, markershape = :circle, aspect_ratio = 2 * Δx/Δy, markercolor = :magenta)
inds = 30:i
p2 = plot(les.t[inds] ./ 86400, sqrt.(loss_default[inds]), legend = :bottomright, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Reference" , color = :green, ylims = (-0.0001, 0.017), xlims = (0.0, 8.0), title = "Error in Time")
# p2 =  plot!(les.t[inds] ./ 86400, sqrt.(loss_optimal[inds]), legend = :topleft, xlabel = "days", ylabel = "Error " * celsius, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, label = "Mode" , color = :magenta, title = "Error in Time")

p = plot(p1, p2)
end
gif(anim, pwd() * "/error_animation_fixed_reference.gif", fps = 15)
mp4(anim, pwd() * "/error_animation_fixed_reference.mp4", fps = 15)

##
# siam_mpe
xup = 19.45
xdown = 19.2
yup = 0.0
ydown = -75
Δy = abs(yup-ydown)
Δx = abs(xup -xdown)
i = length(les.t)
p1 = plot(les.T[:,i], les.z, label = "LES", legend = false, 
ylabel = "depth [m]", xlabel = "Temperature " * celsius, 
xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, 
aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, 
gridalpha = 0.25, framestyle = :box, title = "LES", linewidth =3 )
j = 1
𝑪 = parameter_list[j]
loss_default = ℒᵗ(𝑪)
Tᵖ = 𝒢(𝑪)
daystring = @sprintf("%.1f", les.t[i] ./ 86400)
p2 = scatter(Tᵖ[:,i], zᵖ,  markersize = 6, 
grid = true, gridstyle = :dash, gridalpha = 0.25, 
framestyle = :box, markerstrokealpha = 0.0,  
markerstrokewidth = 0.0, markershape = :square,
 aspect_ratio = 2 * Δx/Δy, markercolor = :green, title = "KPP", 
xlims = (xdown, xup), ylims = (ydown, yup), legend = false,
 ylabel = "depth [m]", xlabel = "Temperature " * celsius)
p = plot(p2, p1)
savefig(p, pwd() * "/side_by_side.pdf")

##
# mcmc example
resolution = resolutions[1]
case = cases[1]
# construct filename
filename = pwd() * "/LES/" * case * "_profiles.jld2"
# construct default loss function
N = resolution[1]
Δt = resolution[2]
ℒ = closure_default_loss_function(filename, N = N, Δt = Δt)
# choose default parameters
optimal_𝑪 = copy(default_𝑪)
resolution_label = "_res_" * string(resolution[1])
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_optima.jld2"
mcmc_data = jldopen(filename, "r")
initial_𝑪 = mcmc_data["parameter"]
ℒ⁰ = mcmc_data["loss"]
close(mcmc_data)

# gif for fun, takes a little while to run
const factor = 1
inverse_factor = false
filename = filename = pwd() * "/mcmc_data/" * "toy_example_" *  string(factor) * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
proposal_chain = mcmc_data["proposal_𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
acceptance_rate = sum(e1 .== e2) / length(e1)
println("the acceptance rate was")
println(acceptance_rate)
indmin = argmin(e1)
close(mcmc_data)
index1 = 3
index2 = 4
bools = e1 .< minimum(e1) * 2
tmp_ind = argmax(bools)
if factor > 1
    tmp_ind = 80
    tmp_ind = 830
else
    tmp_ind = 100
    tmp_ind = 1037
end
bins = 200
Cᴿ = 0.3
chain[4, :] *= Cᴿ
proposal_chain[4, :] *= Cᴿ
right_bounds[4] = 4.0
using LinearAlgebra
tail_ind = 10
final_index = 1200
anim = @animate for i in 1:1:final_index
    p1 = histogram2d(chain[index1, tmp_ind:end], chain[index2, tmp_ind:end], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, legend = false, color = cgrad(:blues, rev = false), xtickfont=font(18), ytickfont=font(18), xguidefontsize=18, yguidefontsize = 18, legendfontsize = 18)
    # Starting value
    scatter!(chain[index1, 1:1], chain[index2, 1:1], shape = :circle, color = :blue, label = "starting value", markersize= 15, opacity = 0.5)
    # Ending Value
    completionstring = @sprintf("%.1f", i /final_index * 100)
    scatter!(initial_𝑪[index1, 1:1], initial_𝑪[index2, 1:1] .* Cᴿ, shape = :star, color = :green, label = "optimal value", legend = :topright, markersize= 25, legendfont = font("Times new roman", 13), opacity = 1.0, title = completionstring * "%  complete")
    # another way to accomplish similar things is with
    # marker_z = (+), color = :bluesreds
    # see http://docs.juliaplots.org/latest/generated/plotly/#plotly-ref35-1
    ω = i / tmp_ind / 8 * 4.5
    p1 = scatter!(chain[index1, i:i+tail_ind], chain[index2, i:i+tail_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 0.5, RGB(0.0, 0.0, 1.0), stroke(0.1, 0.1, :black, :dot)), label = false)

    scatter!(chain[index1, i+1+tail_ind:i+1+tail_ind], chain[index2, i+1+tail_ind:i+1+tail_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 1.0, RGB(0.0, 0.0, 1.0), stroke(0.1, 0.1, :black, :dot)), label = "chain")
    # proposal parameter
    accepted = norm(proposal_chain[:, i+2+tail_ind] - chain[:, i+2+tail_ind]) < eps(1.0)
    if accepted
        color = RGB(0.0, 1.0, 0.0)
    else
        color = RGB(1.0, 0.0, 0.0)
    end
    scatter!(proposal_chain[index1, i+2+tail_ind:i+2+tail_ind], proposal_chain[index2, i+2+tail_ind:i+2+tail_ind], xlabel = parameter_dictionary[index1], ylabel = parameter_dictionary[index2], bins = bins, normalize = true, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (left_bounds[index1], right_bounds[index1]), ylims = (left_bounds[index2], right_bounds[index2]), marker = (:hexagon, 6, 1.0, color, stroke(1, 1.0, :black, :dot)), label = "proposal")

end
gif(anim, pwd() * "/rwmcmc_example.gif", fps = 60)
mp4(anim, pwd() * "/rwmcmc_example.mp4", fps = 60)