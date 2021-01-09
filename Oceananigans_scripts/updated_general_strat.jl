# for version 0.36.0
using Printf, Oceananigans
using Oceananigans.BoundaryConditions,
      Oceananigans.Fields,
      Oceananigans.OutputWriters,
      Oceananigans.Diagnostics,
      Oceananigans.Utils,
      Oceananigans.AbstractOperations
using CUDA
CUDA.allowscalar(true)
# Architecture and float type
arch = GPU()
FT   = Float64

write_output = true 
output_interval     = 600   # [seconds]
checkpoint_interval = 86400 # [seconds]
const scale = 64;
filename_1 = "general_strat_" * string(scale)

# Simulation time
days = 0.5  * scale
end_time = day * days

const Lx = Ly = Lz = 100 # [meters]

# Rough resolution
const Nx = Ny = Nz = 2^7

# Domain
topology = (Periodic, Periodic, Bounded)
grid = RegularCartesianGrid(topology=topology, 
                            size=(Nx, Ny, Nz), 
                            x=(0, Lx), y=(0, Ly), z=(-Lz, 0))

# Coriolis parameter 
const f = -1e-4 # [s⁻¹]
coriolis = FPlane(FT, f=f)
# Constants
const α = 2e-4     # [K⁻¹]
const g = 9.80665  # [m/s²]
const cₚ = 4000.0  # Specific heat capacity of seawater at constant pressure [J/(kg·K)]
const ρ₀ = 1027.0  # Density of seawater [kg/m³]
# Surface Boundary Conditions
bc_params = ()
const τ = -0.0           # Wind stress [N]
const Q = 100.0          # surface flux [W/m²]
# Surface Fluxes
const Φᵘ = τ / ρ₀        # Momentum Flux [m² / s²]
const Qᵇ = Q / (ρ₀ * cₚ) # Buoyancy Flux [m²/s³]

# Initial Stratification
const ∂b∂z = α * g *  0.01 / 16 * scale

# Boundary Condition Implementation
@inline wind_stress(x, y, t, p)  = Φᵘ
@inline surface_flux(x, y, t, p) = Qᵇ
@inline bottom_strat(x, y, t, p) = ∂b∂z

# Boundary Conditions
# Buoyancy
top_b_bc = BoundaryCondition(Flux, surface_flux, parameters = bc_params)
bottom_b_bc = BoundaryCondition(Gradient, bottom_strat, bc_params)
b_bcs = TracerBoundaryConditions(grid, top = top_b_bc, bottom = bottom_b_bc)
# Zonal Velocity
BoundaryCondition(Flux, wind_stress, parameters = bc_params)
u_bcs = UVelocityBoundaryConditions(grid, top = top_u_bc)

# boundary conditions
bcs = (b = b_bcs,  u = u_bcs)

# Initial Conditions
ε(σ) = σ * randn()
# make initial condition same as relaxation to northern wall
B₀(x, y, z) = α * g * 20 + ∂b∂z * z +  α * g * ε(1e-10) * exp(4z/Lz)


# Linear Equation of state
eos = LinearEquationOfState(FT, α=α, β=0)
buoyancy = BuoyancyTracer()

# Closure
closure = AnisotropicMinimumDissipation(FT)

# checkpointing
searchdir(path, key) = filter(x -> occursin(key, x), readdir(path))
checkpoints = searchdir(pwd(), filename_1 * "_checkpoint_iteration")
if length(checkpoints) > 0
    checkpointed = true
    checkpoint_filepath = joinpath(pwd(), checkpoints[end])
    @info "Restoring from checkpoint: $checkpoint_filepath"
    model = restore_from_checkpoint(checkpoint_filepath, boundary_conditions = bcs)
else
    checkpointed = false
	model = IncompressibleModel(
	           architecture = arch,
	             float_type = FT,
	                   grid = grid,
	               coriolis = coriolis,
	               buoyancy = buoyancy,
	                closure = closure,
	                tracers = (:b,),
	    boundary_conditions = bcs
	)
end
if !checkpointed
	set!(model, b=B₀)
end

checkpointer = Checkpointer(model, prefix = filename_1 * "_checkpoint", 
                            time_interval = checkpoint_interval, force = true)
##
# Diagnostics 
# Create output writer that writes vertical profiles to JLD2 output files.
# NaN checker will abort simulation if NaNs are produced.
push!(model.diagnostics, NaNChecker(model; frequency=1000, fields=Dict(:w => model.velocities.w)))

# Define horizontal average diagnostics.
 Up = Average(model.velocities.u;       return_type=Array, dims = (1,2,))
 Vp = Average(model.velocities.v;       return_type=Array, dims = (1,2,))
 Wp = Average(model.velocities.w;       return_type=Array, dims = (1,2,))
 Tp = Average(model.tracers.T;          return_type=Array, dims = (1,2,))
 Sp = Average(model.tracers.S;          return_type=Array, dims = (1,2,))
 νp = Average(model.diffusivities.νₑ;   return_type=Array, dims = (1,2,))
κTp = Average(model.diffusivities.κₑ.T; return_type=Array, dims = (1,2,))
κSp = Average(model.diffusivities.κₑ.S; return_type=Array, dims = (1,2,))

u = model.velocities.u
v = model.velocities.v
w = model.velocities.w
T = model.tracers.T
S = model.tracers.S

uu = Average(u*u, model; return_type=Array, dims = (1,2,))
vv = Average(v*v, model; return_type=Array, dims = (1,2,))
ww = Average(w*w, model; return_type=Array, dims = (1,2,))
uv = Average(u*v, model; return_type=Array, dims = (1,2,))
uw = Average(u*w, model; return_type=Array, dims = (1,2,))
vw = Average(v*w, model; return_type=Array, dims = (1,2,))
wT = Average(w*T, model; return_type=Array, dims = (1,2,))
wS = Average(w*S, model; return_type=Array, dims = (1,2,))
vshear = Average(∂z(u)^2 + ∂z(v)^2, model; return_type=Array, dims = (1,2))
# Create output writer that writes vertical profiles to JLD2 output files.
profiles = Dict(
     :u => model -> Up(model),
     :v => model -> Vp(model),
     :w => model -> Wp(model),
     :T => model -> Tp(model),
     :S => model -> Sp(model),
    :nu => model -> νp(model),
:kappaT => model -> κTp(model),
:kappaS => model -> κSp(model),
    :uu => model -> uu(model),
    :vv => model -> vv(model),
    :ww => model -> ww(model),
    :uv => model -> uv(model),
    :uw => model -> uw(model),
    :vw => model -> vw(model),
    :wT => model -> wT(model),
  :wS => model -> wS(model),
  :vshear => model -> vshear(model)
)

profile_writer = JLD2OutputWriter(model, profiles; 
prefix = "general_strat_" * string(scale) * "_profiles",
time_interval=output_interval, verbose=true)

push!(model.output_writers, profile_writer)
###
if !checkpointed
	Δt= 1.0
else
	Δt = 3.0
end
Δt_wizard = TimeStepWizard(cfl=0.3, Δt = Δt, max_change=1.1, max_Δt= 3.0)
cfl = AdvectiveCFL(Δt_wizard)

# Take Ni "intermediate" time steps at a time before printing a progress
# statement and updating the time step.
Ni = 1000

function print_progress(simulation)
    model = simulation.model
    i, t = model.clock.iteration, model.clock.time

    progress = 100 * (model.clock.time / end_time)

    umax = maximum(abs, model.velocities.u.data.parent)
    vmax = maximum(abs, model.velocities.v.data.parent)
    wmax = maximum(abs, model.velocities.w.data.parent)

    @printf("[%05.2f%%] i: %d, t: %.2e days, umax: (%6.3e, %6.3e, %6.3e) m/s, CFL: %6.4e, next Δt: %.1e s\n",
    	    progress, i, t / day, umax, vmax, wmax, cfl(model), Δt_wizard.Δt)
end

simulation = Simulation(model, Δt=Δt_wizard, stop_time=end_time, progress=print_progress, iteration_interval=Ni)
simulation.output_writers[:checkpoint] = checkpointer
###
run!(simulation)
write_output(model, checkpointer)