using JLD2
files = [
    pwd() * "/three_layer_constant_fluxes_hr48_Qu0.0e+00_Qb1.2e-07_f1.0e-04_Nh256_Nz128_free_convection_averaged_statistics.jld2",
    pwd() * "/three_layer_constant_fluxes_hr48_Qu3.0e-04_Qb1.0e-07_f1.0e-04_Nh256_Nz128_weak_wind_strong_cooling_averaged_statistics.jld2",
    pwd() * "/three_layer_constant_fluxes_hr48_Qu0.0e+00_Qb1.2e-07_f1.0e-04_Nh256_Nz128_free_convection_statistics.jld2"
]
filename = files[3]
file = jldopen(filename)

αg = 0.00196133

keys(file["timeseries"]["u"])
u = []
v = []
b = []
e = []
t = []
for key in keys(file["timeseries"]["u"])
    push!(u,file["timeseries"]["u"][key]*1.0)
    push!(v,file["timeseries"]["v"][key]*1.0)
    push!(b,file["timeseries"]["b"][key]*1.0)
    push!(e,file["timeseries"]["e"][key]*1.0)
    push!(t,file["timeseries"]["t"][key]*1.0)
end
Δt = (t[2] - t[1]) # roughly 10 minute timesteps
zC = file["grid"]["zC"]
H = abs((zC[3] + zC[4])/2)
θ_top = file["boundary_conditions"]["θ_top"] # temperature flux
θ_bottom = file["boundary_conditions"]["θ_bottom"] # temperature gradient
u_top = file["boundary_conditions"]["u_top"]
u_bottom = file["boundary_conditions"]["u_bottom"]
Qb = Meta.parse(split(filename, "_Qb")[2][1:7])
# αg = Qb / θ_top when θ_top isn't zero
b_bottom = θ_bottom * αg
coriolis = Meta.parse(split(filename, "_f")[3][1:7])
##
N = length(b[1])
# Build the model with a Backward Euler timestepper
# KPP.Model, TKEMassFlux
model = KPP.Model(N=N, H=H, stepper=:BackwardEuler)

model.solution.T.data[1:N] .= b[1][:]
model.solution.U.data[1:N] .= u[1][:]
model.solution.V.data[1:N] .= v[1][:]
model.bcs.T.top = FluxBoundaryCondition(Qb)
model.bcs.T.bottom = GradientBoundaryCondition(b_bottom)
model.bcs.U.top = FluxBoundaryCondition(u_top)

i = 13
run_until!(model, Δt, t[i])


plot(model.solution.T.data[1:N], model.grid.zc, label = "KPP" )
plot!(b[i][:], model.grid.zc, legend = :bottomright, label = "LES")

mean(model.solution.T.data[1:N])
##
mean(b[1]) - Qb * t[i] / H # theoretical value
mean(model.solution.T.data[1:N]) # KPP
mean(b[i]) # LES

# should be constant
for i in 2:length(b)
    println(mean(b[i-1]-b[i]) / (t[i] - t[i-1]) * H)
end
gr(size = (400,300))
plot(t, mean.(b), xlabel = "time [s]", ylabel = "average b [m/s²]", legend = false, title = "Mean Buoyancy vs time")
