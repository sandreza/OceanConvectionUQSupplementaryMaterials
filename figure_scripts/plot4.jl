include("../src/LocalOceanUQSupplementaryMaterials.jl")
include("../scripts/utils.jl")
include("../figure_scripts/utils.jl")
# uncertainty propagation figures

display_plot = false
save_figures = true
pyplot()
case = cases[1]
resolution = resolutions[1]

resolution_label = "_res_" * string(resolution[1])
# get posterior data
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_uncertainty_propagation.jld2"
histogram_data = jldopen(filename, "r")
h1 = histogram_data["h1"]
close(histogram_data)
# get prior data
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_prior_uncertainty_propagation.jld2"
histogram_data = jldopen(filename, "r")
h1p = histogram_data["h1"]
close(histogram_data)
# get domain data
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_domain.jld2"
domain_data = jldopen(filename, "r")
zᵖ = domain_data["zᵖ"]
ϕrange = domain_data["ϕrange"]
ϕrange = (ϕrange[1:(end-1)] + ϕrange[2:end] )./2
close(domain_data)

max_ind = length(h1)

mat = uq_prop_mat(h1, index = max_ind-20)
cmap = heatmap(ϕrange, zᵖ, log.(mat')/log(10),  color =:fire, background_color= :white, xlabel = "Temperature " * celsius, clims = (-5,-1.0), title="Posterior Ensemble", ylims = (-80,0), xlims = (19.2,19.5), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, yaxis = false, left_margin = 10mm);
if display_plot
    display(cmap)
end

matp = uq_prop_mat(h1p, index = max_ind-20)
cmap_prior = heatmap(ϕrange, zᵖ, log.(matp')/log(10),  color =:fire, background_color= :white, xlabel = "Temperature " * celsius, ylabel = "Depth [m]", clims = (-5,-1.0), title= "Prior Ensemble", ylims = (-80,0), xlims = (19.2,19.5), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, colorbar = false);
if display_plot
    display(cmap_prior)
end

p1 = plot(cmap_prior, cmap);
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_4.pdf")
end

###
# loop for gif purposes
anim = @animate for i in 1:1:max_ind
    mat = uq_prop_mat(h1, index = i)
    matp = uq_prop_mat(h1p, index = i)

    cmap_posterior = heatmap(ϕrange, zᵖ, log.(mat')/log(10),  color =:fire, background_color= :white, xlabel = "Temperature " * celsius, ylabel = "Depth [meters]", clims = (-5,-0.5), title= "Posterior Ensemble", ylims = (-80,0))

    cmap_prior = heatmap(ϕrange, zᵖ, log.(matp')/log(10),  color =:fire, background_color= :white, xlabel = "Temperature " * celsius, ylabel = "Depth [meters]", clims = (-5,-0.5), title= "Prior Ensemble", ylims = (-80,0))
    if display_plot
        display(plot(cmap_prior,cmap_posterior))
    end
end

if save_figures == true
    gif(anim, pwd() * "/figures/figure_4_animation.gif", fps = 15)
end
