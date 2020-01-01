include("./src/LocalOceanUQSupplementaryMaterials.jl")
include("./scripts/utils.jl")
include("./figure_scripts/utils.jl")

# compromise functions

save_figures = true

resolution = resolutions[1]
# define things for forward map
N = resolution[1]
Δt = resolution[2]
zᵖ = zeros(N)

ind_case_1 = 1
ind_case_2 = 2

# get LES 1
case = cases[ind_case_1]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les1 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les1.t)
𝒢1 = CoreFunctionality.closure_free_convection_flexible(N, Δt, les1, subsample = subsample, grid = zᵖ, power = 1.0)
ℒ1 = CoreFunctionality.closure_T_nll(𝒢1, les1; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )
filename = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
indmin1 = argmin(e1)
chain1 = mcmc_data["𝑪"]
𝑪1 = chain1[:,indmin1]
𝑪1 = median(chain1, dims = 2)[:]
println("median 1 is $𝑪1")
### now fo rcase 2
case = cases[ind_case_2]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
les2 = CoreFunctionality.OceananigansData(filename)
subsample = 1:1:length(les2.t)
𝒢2 = CoreFunctionality.closure_free_convection_flexible(N, Δt, les2, subsample = subsample, grid = zᵖ, power = 1.0)
ℒ2 = CoreFunctionality.closure_T_nll(𝒢2, les2; weight = 1, subsample = subsample, series=true, power = 2, f1 = mean, f2 = maximum )
filename = pwd() * "/mcmc_data/" * case * resolution_label  * "_flexible_new" * "_mcmc.jld2"
mcmc_data = jldopen(filename, "r")
chain = mcmc_data["𝑪"]
e1 = mcmc_data["ε"]
e2 = mcmc_data["proposal_ε"]
indmin2 = argmin(e1)
chain2 = mcmc_data["𝑪"]
𝑪2 = chain2[:,indmin2]
𝑪2 = median(chain2, dims = 2)[:]
println("median 2 is $𝑪2")
# create other compromise distribution
chain4 = combine(chain1, chain2)
𝑪4 = median(chain4,dims=2)[:] #other choice of compromise

NN1 = 𝑪1[6]
NN2 = 𝑪2[6]

###
#=
case_range1 = 1:(10^3-1)
case_range2 = copy(case_range1)
case_range4 = (5 * 10^2 +2) : (1 * 10^3 - 1 + 5 * 10^2)
index = 5
histogram(chain1[index,case_range1], normalize = true, alpha = 0.4, xlims = (left_bounds[index], right_bounds[index]), xlabel = parameter_dictionary[index], bins = 50, legend = true, ylabel = "pdf", label = "Case 1")
histogram!(chain2[index,case_range2], normalize = true, alpha = 0.4, label = "Case 2")
histogram!(chain4[index,case_range4], normalize = true, alpha = 0.4, label = "Compromise")
###
p = marginal_pdfs(chain1[:,case_range1], chain2[:,case_range2], left_bounds, right_bounds, parameter_dictionary)
plot(p...)
p1 = plot(p[4])

if save_figures == true
    savefig(p1, pwd() * "/figures/figure_9_distributions.png")
end
=#
###
# parameters to loop over
# show everything in case 1 scenario
p_case1 = []
labels = ["Median 1", "Median 2",  "Compromise"]
#𝑪1[1] = 1e-4
#𝑪2[2] = 3.5
#𝑪2[5] = 0.375
#𝑪2[3] = 10
parameter_list =[𝑪1, 𝑪2, 𝑪4]
# the stratification changes
𝑪1[6] = NN1
𝑪2[6] = NN1
𝑪4[6] = NN1
maxT = maximum(les1.T[:,end] .+ 0.1)
minT = minimum(les1.T[:,end])
plot()
for j in 1:3
    𝑪 = parameter_list[j]
    loss = maximum( ℒ1(𝑪) )
    Tᵖ = 𝒢1(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les1.T[:,end], les1.z, label = "LES 1", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature [C]", ylims = (-100, 0), xlims = (minT, maxT))
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 1 = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case1,p1)
end
p1 = plot(p_case1...)
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_9_case_1.png")
end
###
# show everything in case 2 scenario
p_case2 = []
labels = ["Median 1", "Median 2", "Compromise"]
parameter_list =[𝑪1, 𝑪2, 𝑪4]
𝑪1[6] = NN2
𝑪2[6] = NN2
𝑪4[6] = NN2
plot()
maxT = maximum(les2.T[:,end] .+ 0.1)
minT = minimum(les2.T[:,end])
xlims = (minT,maxT)
for j in 1:3
    𝑪 = parameter_list[j]
    loss = maximum(ℒ2(𝑪))
    Tᵖ = 𝒢2(𝑪)
    loss_string = @sprintf("%.1e", sqrt(loss))
    p1 = plot(les2.T[:,end], les2.z, label = "LES 2", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature [C]", ylims = (-100,0), xlims = (minT,maxT))
    p1 = scatter!(Tᵖ[:,end], zᵖ, label = labels[j], title = "Error 2 = " * loss_string * " [C]", grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)
    display(p1)
    push!(p_case2,p1)
end
p2 = plot(p_case2...)
if save_figures == true
    savefig(p2, pwd() * "/figures/figure_11_case_2.png")
end
###
# Show median values for optimal in both cases
p1 = plot(p_case1[1], p_case2[2])
if save_figures == true
    savefig(p1, pwd() * "/figures/figure_11_alternative.png")
end

# Show median values in opposite cases
p2 = plot(p_case1[2], p_case2[1])
if save_figures == true
    savefig(p2, pwd() * "/figures/figure_11.png")
end

# Show them together
p4 = plot(p_case1[1], p_case2[2], p_case1[2], p_case2[1])
if save_figures == true
    savefig(p4, pwd() * "/figures/figure_11_alternative3.png")
end

# Show compromise values in both cases
p5 = plot(p_case1[3], p_case2[3])
if save_figures == true
    savefig(p5, pwd() * "/figures/figure_11_alternative4.png")
end

# special = p3
#plot(special, p4)
#plot(tmp)


####
using OceanTurb
# Build the model with a Backward Euler timestepper
𝑪 = [0.11803164331592443, 3.7246545857676954, 0.35191154207167974, 6.225750233165317]
parameters = KPP.Parameters( CSL = 𝑪[1], CNL = 𝑪[2], Cb_T = 𝑪[3], CKE = 𝑪[4])
constants = Constants(Float64; α = les.α , β = les.β, ρ₀= les.ρ, cP=les.cᵖ, f=les.f⁰, g=les.g)
N = 16
Δt = 600.0
model = KPP.Model(N=N, L=les.L, stepper=:BackwardEuler, constants = constants) #, parameters = parameters)

# get average of initial condition of LES
T⁰ = avg(les.T⁰, N)
# set equal to initial condition of parameterization
model.solution.T[1:N] = copy(T⁰)
# Set boundary conditions
model.bcs.T.top = FluxBoundaryCondition(les.top_T)
model.bcs.T.bottom = GradientBoundaryCondition(les.bottom_T)
# set aside memory
if subsample != 1
    time_index = subsample
else
    time_index = 1:length(les.t)
end
Nt = length(les.t[time_index])
𝒢 = zeros(Nt)
𝒢1 = zeros(N, Nt)

# loop the model
ti = collect(time_index)
for i in 1:Nt
    t = les.t[ti[i]]
    run_until!(model, Δt, t)
    @. 𝒢1[:,i] = model.solution.T[1:N]
    𝒢[i] = model.state.h
end
# 60 - 80
z = collect(model.grid.zc)
scatter(ti, 𝒢)
scatter(𝒢1[:, 80], z, legend = false)



###
calc_prior = false
const number_of_ensembles = 10^6
const skip = 100
for case in cases[1:1]
    for resolution in resolutions[1:1]
        # construct filename
        filename = pwd() * "/LES/" * case * "_profiles.jld2"
        # load les
        les = CoreFunctionality.OceananigansData(filename)
        # construct default loss function
        N = resolution[1]
        Δt = resolution[2]
        # define the forward map
        zᵖ = zeros(N)
        #calculate every hour
        subsample = 1:6:length(les.t)
        # define the forward map
        𝒢 = CoreFunctionality.closure_free_convection_ml_depth(N, Δt, les, subsample = subsample, grid = zᵖ)
        println("-------------------")
        println("For case $case ")
        println("and resolution " * string(resolution[1]))
        resolution_label = "_res_" * string(resolution[1])
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_mcmc.jld2"
        mcmc_data = jldopen(filename, "r")
        𝑪 = mcmc_data["𝑪"]
        close(mcmc_data)
        Φ = 𝒢(𝑪)
        ϕmin = 0.0
        ϕmax = 100.0
        Δϕ = (ϕmax - ϕmin) / 1000
        ϕrange = collect(ϕmin:Δϕ:ϕmax)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_uncertainty_propagation_ml.jld2"
        CoreFunctionality.propagate_uncertainty(𝑪[:,1:skip:number_of_ensembles], 𝒢, field_range = ϕrange, filename = filename)
        filename = pwd() * "/mcmc_data/" * case * resolution_label * "_domain.jld2"
        println("done with posterior")
        if calc_prior
            # Now do it for the prior distribution
            filename = pwd() * "/mcmc_data/" * "prior" * "_mcmc.jld2"
            mcmc_data = jldopen(filename, "r")
            𝑪 = mcmc_data["𝑪"]
            close(mcmc_data)
            filename = pwd() * "/mcmc_data/" * case * resolution_label * "_prior" * "_uncertainty_propagation_ml.jld2"
            CoreFunctionality.propagate_uncertainty(𝑪[:,1:skip:number_of_ensembles], 𝒢, field_range = ϕrange, filename = filename)
            println("done with prior")
        end
    end
end

###
save_figures = false

case = cases[1]
resolution = resolutions[1]

resolution_label = "_res_" * string(resolution[1])
# get posterior data
# filename = pwd() * "/mcmc_data/" * case * resolution_label * "_prior_uncertainty_propagation_ml.jld2"
filename = pwd() * "/mcmc_data/" * case * resolution_label * "_uncertainty_propagation_ml.jld2"
histogram_data = jldopen(filename, "r")
h1 = histogram_data["h1"]
close(histogram_data)

m = length(h1[100][1].weights)
n = length(h1)
mat = randn(m,n)
for i in 1:1:193
    ind = i
    normalization_constant = sum(h1[ind][1].weights)
    tmp = h1[ind][1].weights / normalization_constant
    @. mat[:,i] = tmp
end
ϕmin = 0.0
ϕmax = 100.0
Δϕ = (ϕmax - ϕmin) / 1000
ϕrange = collect(ϕmin:Δϕ:ϕmax)
ϕrange = (ϕrange[2:end] + ϕrange[1:end-1])./2
trange = les.t[1:6:length(les.t)] ./ 86400
p1 = heatmap(trange, ϕrange, log.(mat .+ eps(1.0))./log(10), ylabel = "Mixed Layer Depth [m]", xlabel = "days", title = "Mixed Layer Depth Uncertainty", clims = (-3, -0), ylims = (0, 75))
display(p1)
savefig(p1, pwd() * "/figures/ml_figure.png")
println("at time ")
println(trange[end-24])
println("days")
bools = h1[end-30][1].weights .> 8

seen = ϕrange[bools]
println(seen)
max_seen = maximum(seen[1:end])
min_seen = minimum(seen[1:end])
println("The largest mixed layer depth is $max_seen")
println("The smallest mixed layer depth is $min_seen")

timelabel = @sprintf("%.1f", trange[end])
tmp = ones(length(h1[end][1].weights))
p1 = histogram(ϕrange, weights = h1[end][1].weights, bins = 200, norm = true, xlims = (65, 75), ylabel = "probability", xlabel = "mixed layer depth", title = "Mixed layer depth uncertainty at " * timelabel * " days", legend = false)
savefig(p1, pwd() * "/figures/ml_figure_alternative.png")
