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
# plot LES data and Initial Condition
xup = 20.01
xdown = 19.2
yup = 0.0
ydown = -75
Δy = abs(yup-ydown)
Δx = abs(xup -xdown)

case = cases[1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)


p1 = plot(20 .+ les.z .* 0.01, les.z, label = "time 0",
ylabel = "depth [m]", xlabel = "Temperature " * celsius, 
xlims = (xdown, xup), ylims = (ydown, yup), color = :green, 
aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, 
gridalpha = 0.25, framestyle = :box, linewidth = 3)
p1 = plot!(les.T[:,end], les.z, label = "time t", legend = :bottomright, 
ylabel = "depth [m]", xlabel = "Temperature " * celsius, 
xlims = (xdown, xup), ylims = (ydown, yup), color = :blue, 
aspect_ratio = 2 * Δx/Δy, grid = true, gridstyle = :dash, 
gridalpha = 0.25, framestyle = :box, linewidth = 3)
display(p1)
savefig(p1, pwd() * "/initial_and_final.pdf")
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

##
# Mixed layer depth animations

using JLD2
include(pwd() * "/src/LocalOceanUQSupplementaryMaterials.jl")
include(pwd() * "/scripts/utils.jl")
include(pwd() * "/figure_scripts/utils.jl")

# les analysis
# derivative function
function δ(z, Φ)
    m = length(Φ)-1
    Φz = ones(m)
    for i in 1:m
        Φz[i] = (Φ[i+1]-Φ[i])/(z[i+1]-z[i])
    end
    return Φz
end



animation = false
#case = cases[18]
case = "dns_old"
# case = "dns2"
case = "rdns"
case = "rdns_2"
case = "dns"
case = cases[1]
newcases = ["dns_old", "dns_2", "ndns", "rdns", "rdns_2"]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)

# get vertical shear
les_data = jldopen(filename, "r")
les_keys = keys(les_data)
timeseries_keys = keys(les_data["timeseries"]["t"])
list = keys(les_data["timeseries"])
timeseries_keys = keys(les_data["timeseries"]["t"])
# size of arrays
Nz = length(collect(les_data["grid"]["zC"]))
Nt = length(timeseries_keys)
vshear = zeros(Nz,Nt)
# get vertical shear
if case in cases2
    println("extracting vertical shear")
    for j in 1:Nt
        # Fields
        key = timeseries_keys[j]
        @. vshear[:,j] = les_data["timeseries"]["vshear"][key][2:(end-1)]
    end
elseif case in newcases
    println("extracting vertical shear")
    for j in 1:Nt
        # Fields
        key = timeseries_keys[j]
        @. vshear[:,j] = les_data["timeseries"]["vshear"][key][2:(end-1)]
    end
else
    nothing
end
close(les_data)

Qᵇ = les.α * les.g * les.top_T
N² = les.α * les.g * les.bottom_T
Nt = length(les.t)

maxww = maximum(les.ww, dims =1)[:]
h = randn(Nt)
h1 = rand(Nt)
h2 = randn(Nt)
ΔB = randn(Nt)
ΔB2 = randn(Nt)
Vᵗ = randn(Nt)
we = randn(Nt)
ww_top = randn(Nt)
ww_base = randn(Nt)
ww_star = randn(Nt)
Nᵉ = randn(Nt)
TKE_base = randn(Nt)
ΔTKE = randn(Nt)
𝒮 = randn(Nt)
for i in 1:Nt
    B = les.α * les.g * les.T[:,i]
    Bz = δ(les.z, B)
    mBz = maximum(Bz)
    hⁱ = argmin(les.wT[:,i])
    hⁱ = argmax(Bz)
    h[i] = -les.z[hⁱ]
    regind = maximum([1, hⁱ-10])
    ΔB[i] = mean(B[(end-10):end]) - B[hⁱ]
    Vᵗ[i] = ( h[i] * ΔB[i] ) / ( h[i] * sqrt(N²)  * (h[i] * Qᵇ)^(1/3) )
    we[i] = -(les.α * les.g * les.wT[hⁱ,i]) ./ ΔB[i]

    Nᵉ[i] = sqrt(mBz)
    tt = (2*N² + mBz)/3
    bools = Bz .> tt
    bools2 = Bz .> (N² + 0)/4
    ind = collect(1:length(Bz))
    cand = ind[bools]
    zA = (les.z[1:(end-1)] + les.z[2:end] )./2
    bools3 = zeros(Bool,length(bools))
    @. bools3[hⁱ:end] = bools[hⁱ:end]
    h2[i] = -minimum(zA[bools])
    h1[i] = -maximum(zA[bools2])
    h1ⁱ = argmax(zA[bools])
    h2ⁱ = argmin(zA[bools])
    TKE_base[i] = les.uu[hⁱ,i] + les.vv[hⁱ,i] + les.ww[hⁱ,i]
    TKE_base[i] = maximum(les.uu[:,i] + les.vv[:,i] + les.ww[:,i])
    # TKE_base[i] = les.uu[h1ⁱ,i] + les.vv[h1ⁱ,i] + les.ww[h1ⁱ,i]
    ΔTKE[i] = les.uu[h1ⁱ,i] + les.vv[h1ⁱ,i] + les.ww[h1ⁱ,i] - les.uu[h2ⁱ,i] + les.vv[h2ⁱ,i] + les.ww[h2ⁱ,i]
    𝒮[i] = maximum(vshear[1:h1ⁱ, i])
    ww_top[i] = les.ww[h1ⁱ, i]
    ww_base[i] = les.ww[h2ⁱ, i]
    ww_star[i] = (Qᵇ * h1[i])^(2/3)
    ΔB2[i] = mean(B[(end-10):end]) - B[h2ⁱ]
end
field = les.T
max_field = maximum(field)
min_field = minimum(field)
animation = false
end_ind = floor(Int, Nt/1)
begin_ind = 30
anim = @animate for i in begin_ind:10:end_ind
    ϕ = field[:,i]
    time_string = @sprintf("%.1i", les.t[i] ./ 86400)
    p1 = plot(ϕ, les.z, legend = :bottomright, xlims = (min_field, max_field), title = "Day = " * time_string, label = "LES", color = :blue, linewidth = 3)
    #plot!(max_field .+ (les.z .* (max_field - min_field)/les.L), -h[i] .+ (les.z .* 0), label = "h_e", linewidth = 2 )
    #plot!(max_field .+ (les.z .* (max_field - min_field)/les.L), -h1[i] .+ (les.z .* 0), label = "h_1", linewidth = 2 )
    plot!(max_field .+ (les.z .* (max_field - min_field)/les.L), -h2[i] .+ (les.z .* 0), label = "h", linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box , ylabel = "Depth [meters]", xlabel = "Temperature [Celcius]", color = :red)
    e_field = @. les.α * les.g * les.wT[:,i] / Qᵇ * (max_field - min_field) ./ 2 + (min_field + max_field )/2
    #plot!( e_field, les.z,  label = "bouyancy flux")
    #plot!( 0 .* les.z .+ (min_field + max_field )/2 , les.z,  label = "zero bouyancy flux", legend = false)
    t = les.t[begin_ind:30:end_ind]
    p2 = scatter(t ./ 86400, sqrt.(t .* Qᵇ / N² * 3.02), marker = :square, color = :green, label = "Fit", legend = :bottomright, xlabel = "Day", ylabel = "h [meters]", title = "Growth in Time")
    p2 = plot!(les.t[begin_ind:i] ./ 86400 , h2[begin_ind:i], label = "LES", color  = :red, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, linewidth = 3)
    p = plot(p1,p2)
    display(plot(p1, p2))
end
gif(anim, pwd() * "/depth_growth_example.gif", fps = 60)
mp4(anim, pwd() * "/depth_growth_example.mp4", fps = 30)