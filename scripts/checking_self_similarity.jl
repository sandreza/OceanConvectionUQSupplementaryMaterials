include(pwd() * "/scripts/utils.jl")
include(pwd() * "/src/LocalOceanUQSupplementaryMaterials.jl")
include(pwd() * "/src/utils.jl")

##
case = cases[3]
filename = pwd() * "/LES/" * case * "_profiles.jld2"
# load CoreFunctionality and Utils
les = CoreFunctionality.OceananigansData(filename)

i = 1
tmp = avg( les.T[:, i],    16 )
z = avg(les.z, 16)
tmplist = [100 * mean(les.T[:,i] - les.T[:, i + 1]) / (les.t[i+1] - les.t[i]) for i in 1:length(les.t)-1]
# dimensional parameters
# non-dimensionalization
##
α = les.α
g = 9.806
N² = les.bottom_T * α * g
ρ⁰ = les.ρ
cᵖ = 4000
Q = 100
κ = 1e-5
Qᶿ = (Q * α * g) / (cᵖ * ρ⁰ * κ * N²) # Q/cᵖ * ρ⁰ should basically be tmplist  
L = 100
τ = L^2 / κ * (Qᶿ + 1)^(-1)
Θ = (N² * L) / (α * g)
##
base_t = les.t ./ τ
plot(tmp ./ Θ,z ./ L)

base_tmp = copy((les.T .- les.T[1]) ./ Θ)

##
anim = @animate for (i,times) in enumerate(base_t)
    println(i)
    println(times)
    # tmptmp = avg(base_tmp[:, i], 16)
    tmptmp = base_tmp[:, i]
    p1 = plot( tmptmp, les.z ./ L, ylims = (-1, 0), xlims = (0, 1), label = "general_strat_1", legend  = :bottomright)

for case in cases[4:3:end]
    filename = pwd() * "/LES/" * case * "_profiles.jld2"
    # load CoreFunctionality and Utils
    les = CoreFunctionality.OceananigansData(filename)
    z = avg(les.z, 16)
    # dimensional parameters
    # non-dimensionalization
    ##
    α = les.α
    g = 9.806
    N² = les.bottom_T * α * g
    ρ⁰ = les.ρ
    cᵖ = 4000
    Q = 100
    κ = 1e-5
    Qᶿ = (Q * α * g) / (cᵖ * ρ⁰ * κ * N²) # Q/cᵖ * ρ⁰ should basically be tmplist  
    L = 100
    τ = L^2 / κ * (Qᶿ + 1)^(-1)
    Θ = (N² * L) / (α * g)
    new_t = les.t ./ τ
    time_index = argmin( abs.(new_t .- times))
    #tmp = avg( les.T[:, time_index] .- les.T[1, 1],    16 )
    #plot!(tmp ./ Θ, z ./ L, label = case)
    tmp = les.T[:, time_index] .- les.T[1, 1]
    plot!(tmp ./ Θ, les.z ./ L, label = case)
end

display(p1)

end
gif(anim, pwd() * "/self_similarity.gif", fps = 60)
