
using OceanTurb

const γ1 = 0.01/2 # related to stratification
const β1  = 0.1 # half the jump in temperature
const ϵ1 = 100.0  # half the length scale over which temperature changes rapidly
const ϵ2 = 400

#θ₀(z) = z * γ1  + β1 * tanh((z+ϵ2) / ϵ1) + 10
θ₀(z) = z * γ1 + 10

wd = pwd()
filename = wd * "/LES/general_strat_1_profiles.jld2"
les = CoreFunctionality.OceananigansData(filename)
Qᵇ = ( 200 / (les.ρ * les.cᵖ) ) * les.α * les.g
N² = γ1 * les.α * les.g
println("The stratification is $(N²)")
N = 64
final_time = 86400 * 30.5 * 6
record_interval = 86400
L = 1000.0
Δt = 100 * 60
# set parameters
𝑪 = [0.10352574433216093, 4.3556149018959, 1.2080069234475657, 7.359241343350497, 0.0 ,0.0]
𝑪 = [0.1, 6.33, 1.36, 7.00, 0.0, 0.0]
parameters = KPP.Parameters( CSL = 𝑪[1], CNL = 𝑪[2], Cb_T = 𝑪[3], CKE = 𝑪[4], CKE2 = 𝑪[5], CKE3 = 𝑪[6])
# Build the model with a Backward Euler timestepper
constants = Constants(Float64; α = les.α , β = les.β, ρ₀= les.ρ, cP=les.cᵖ, f=les.f⁰, g=les.g)
model = KPP.Model(N=N, L=L, stepper=:BackwardEuler, constants = constants, parameters = parameters)
    # Get grid if necessary
zp = collect(model.grid.zc)

# get average of initial condition of LES
T⁰ = θ₀.(zp)
# set equal to initial condition of parameterization
model.solution.T[1:N] = copy(T⁰)
# Set boundary conditions
model.bcs.T.top = FluxBoundaryCondition(Qᵇ / (les.α * les.g))
model.bcs.T.bottom = GradientBoundaryCondition(N²/(les.α * les.g))
# set aside memory
t = collect(0:record_interval:final_time)
Nt = length(t)
𝒢 = zeros(N, Nt)
h1 = zeros(Nt)
# loop the model
ti = collect(time_index)
for i in 1:Nt
    run_until!(model, Δt, t[i])
    @. 𝒢[:,i] = model.solution.T[1:N]
    h1[i] = model.state.h
end

p1 = plot(𝒢[:,1], zp)
println("mixed layer depth 1 is $(h1[end])")
##
𝑪 = [0.132375367537476, 4.207959580871827, 1.4344368071261107, 3.543929859372545, 0.0 ,0.0]
𝑪 = [0.1, 6.33, 1.36, 3.00, 0.0, 0.0]
parameters = KPP.Parameters( CSL = 𝑪[1], CNL = 𝑪[2], Cb_T = 𝑪[3], CKE = 𝑪[4], CKE2 = 𝑪[5], CKE3 = 𝑪[6])
# Build the model with a Backward Euler timestepper

constants = Constants(Float64; α = les.α , β = les.β, ρ₀= les.ρ, cP=les.cᵖ, f=les.f⁰, g=les.g)
model = KPP.Model(N=N, L=L, stepper=:BackwardEuler, constants = constants, parameters = parameters)
    # Get grid if necessary
zp = collect(model.grid.zc)

# get average of initial condition of LES
T⁰ = θ₀.(zp)
# set equal to initial condition of parameterization
model.solution.T[1:N] = copy(T⁰)
# Set boundary conditions
model.bcs.T.top = FluxBoundaryCondition(Qᵇ / (les.α * les.g))
model.bcs.T.bottom = GradientBoundaryCondition(N²/(les.α * les.g))
# set aside memory
t = collect(0:record_interval:final_time)
Nt = length(t)
𝒢2 = zeros(N, Nt)
h2 = zeros(Nt)
# loop the model
ti = collect(time_index)
for i in 1:Nt
    run_until!(model, Δt, t[i])
    @. 𝒢2[:,i] = model.solution.T[1:N]
    h2[i] = model.state.h
end

p2 = plot(𝒢2[:,end], zp)
println("mixed layer depth 2 is $(h2[end])")

plot(p1,p2)


###
plot()
Tmax = maximum(𝒢[:,1])
anim = @animate for i in 1:10:Nt
    p = []
    h1_string = @sprintf("%.1f", h1[i])
    h2_string = @sprintf("%.1f", h2[i])
    p1 = plot(𝒢[:,i], zp, label = "KPP 1", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius)
    plot!(Tmax .+ (zp .* γ1), -h1[i] .+ (zp .* 0), label = "h = " * h1_string * " [m]" , legend = :bottomright, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box )
    push!(p,p1)
    p1 = plot(𝒢2[:,i], zp, label = "KPP 2", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius)
    plot!(Tmax .+ (zp .* γ1), -h2[i] .+ (zp .* 0), label = "h = " * h2_string * " [m]" , legend = :bottomright, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box )
    push!(p,p1)

    display(plot(p...))
end
if save_figures == true
    gif(anim, pwd() * "/figures/figure_10_dynamic.gif", fps = 15)
end

###
p = []
i = Nt
h1_string = @sprintf("%.1f", h1[i])
h2_string = @sprintf("%.1f", h2[i])
p1 = scatter(𝒢[:,i], zp, label = "KPP 1", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius)
plot!(Tmax .+ (zp .* γ1), -h1[i] .+ (zp .* 0), label = "h = " * h1_string * " [m]" , legend = :topleft, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlims = (5, 7.5 ))
p2 = scatter(𝒢2[:,i], zp, label = "KPP 2", legend = :topleft, ylabel = "depth [m]", xlabel = "Temperature " * celsius)
plot!(Tmax .+ (zp .* γ1), -h2[i] .+ (zp .* 0), label = "h = " * h2_string * " [m]" , legend = :topleft, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box , xlims = (5, 7.5) )
p3 = plot(p1,p2)
display(p3)
savefig(p3, pwd() * "/figures/figure_10.png")
