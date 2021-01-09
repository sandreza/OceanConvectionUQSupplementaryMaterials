include(pwd() * "/tests/ensemble_kalman_process.jl")
using LinearAlgebra

rng_seed = 41
Random.seed!(rng_seed)

# Number of synthetic observations from G(u)
n_obs = 1
# Defining the observation noise level
noise_level =  1e-8   
# Independent noise for synthetic observations       
Γy = noise_level * Matrix(I, n_obs, n_obs) 
noise = MvNormal(zeros(n_obs), Γy)

# Loss Function (unique minimum)
function G(u)
    return [sqrt((u[1]-1)^2 + (u[2]+1)^2)]
end

# Loss Function Minimum
u_star = [1.0, -1.0]
y_obs  = G(u_star) + 0 * rand(noise) 

# Define Prior
prior_distns = [Parameterized(Normal(0., sqrt(1))),
                Parameterized(Normal(-0., sqrt(1)))]
constraints = [[no_constraint()], [no_constraint()]]
prior_names = ["u1", "u2"]
prior = ParameterDistribution(prior_distns, constraints, prior_names)
prior_mean = reshape(get_mean(prior),:)
prior_cov = get_cov(prior)

# Calibrate
N_ens = 50  # number of ensemble members
N_iter = 20 # number of EKI iterations
initial_ensemble = construct_initial_ensemble(prior, N_ens;
                                                rng_seed=rng_seed)

ekiobj = EnsembleKalmanProcess(initial_ensemble,
                    y_obs, Γy, Inversion())
##
g_ens = zeros(size(ekiobj.u[1])[1], size(G(ekiobj.u[1][1,:]))[1])


tic = time()
run_eki!(ekiobj, N_iter, G, g_ens, N_ens)
toc = time()
println(tic-tic)
##
for i in eachindex(ekiobj.u)
    p = plot(ekiobj.u[i][:,1], ekiobj.u[i][:,2], seriestype=:scatter, xlims = extrema(ekiobj.u[1][:,1]), ylims = extrema(ekiobj.u[1][:,2]))
    plot!([u_star[1]], xaxis="u1", yaxis="u2", seriestype="vline",
        linestyle=:dash, linecolor=:red, label = false,
        title = "EKI iteration = " * string(i)
        )
    plot!([u_star[2]], seriestype="hline", linestyle=:dash, linecolor=:red, label = "optimum")
    display(p)
    sleep(0.1)
end