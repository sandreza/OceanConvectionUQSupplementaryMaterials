using JLD2
include("./src/LocalOceanUQSupplementaryMaterials.jl")
include("./scripts/utils.jl")
include("./figure_scripts/utils.jl")

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
field = les.ww
max_field = maximum(field)
min_field = minimum(field)
animation = true
end_ind = floor(Int, Nt/1)
if animation
    for i in 1:10:end_ind
        ϕ = field[:,i]
        time_string = @sprintf("%.1i", les.t[i] ./ 86400)
        p1 = plot(ϕ, les.z, legend = :bottomright, xlims = (min_field, max_field), title = "Day = " * time_string, label = "Temperature")
        plot!(max_field .+ (les.z .* (max_field - min_field)/les.L), -h[i] .+ (les.z .* 0), label = "h_e", linewidth = 2 )
        plot!(max_field .+ (les.z .* (max_field - min_field)/les.L), -h1[i] .+ (les.z .* 0), label = "h_1", linewidth = 2 )
        plot!(max_field .+ (les.z .* (max_field - min_field)/les.L), -h2[i] .+ (les.z .* 0), label = "h_2", linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box , ylabel = "Depth [meters]", xlabel = "Temperature [Celcius]")
        e_field = @. les.α * les.g * les.wT[:,i] / Qᵇ * (max_field - min_field) ./ 2 + (min_field + max_field )/2
        #plot!( e_field, les.z,  label = "bouyancy flux")
        #plot!( 0 .* les.z .+ (min_field + max_field )/2 , les.z,  label = "zero bouyancy flux", legend = false)
        plot(p1)
        display(p1)
        println(maximum(ϕ))
    end
end
plot(les.t, maxww)
plot(les.t, h)
plot(les.t ./ 86400, maxww ./ ( (Qᵇ  .* h ).^(2/3) ), ylims = (0, 2))

Ri = (h .* ΔB ) ./ maxww
plot(les.t, h)

plot(les.t[10:end] ./ 86400, Ri[10:end], ylims = (0, 100), title = "Richardson Number")
histogram(Ri[10:end])
median(Ri[10:end])

#plot(les.t ./ 86400, maxww ./ ( N²  .* ( h ).^(2 ) ), ylims = (0,0.01))
plot(les.t, Vᵗ, ylims = (0.5, 2.5))
plot(les.t, ΔB)
histogram(Vᵗ)
plot(les.t, h .* sqrt(N²))
plot!(les.t,(Qᵇ  .* h ).^(1/3) )
scatter(h .* sqrt(N²),  (Qᵇ  .* h ).^(1/3), legend = false)
mean(Vᵗ)
#
plot(les.t[10:end], we[10:end], ylims = (-0.0,0.0002))

scatter(1 ./ Ri[10:end], abs.(we[10:end]) ./ sqrt.(maxww[10:end]), ylims = (0,0.1), xlims = (0,0.1))

plot(les.t, h .* sqrt(N²) ./ (sqrt.(maxww)), ylims = (0, 20))
tmp = Nᵉ ./ sqrt(N²)
scatter(Ri[10:end], tmp[10:end], legend = false)

Δh = h2 - h1

ff(x) = 0.21 + 1.3 / x
#ff(x) = 0.06 + 100 / x^2
scatter(Ri[20:end], Δh[20:end] ./ h1[20:end], ylims = (0,0.3) )
plot!(Ri[20:end], ff.(Ri[20:end]), color = :red, linewidth = 4)

###
scatter(Ri[20:end], Δh[20:end] ./ h1[20:end] .* (Nᵉ[20:end].^2) ./ N², ylims = (0,1.0), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlabel = "Richardson Number", ylabel = "(h2-h1) / h1  * (N_e / N_b)^2" , legend = false)


###
scatter(Ri[20:end], ΔB[20:end] ./ h[20:end] ./ N², ylims = (0,1.3), grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlabel = "Richardson Number", ylabel = "(h2-h1) / h1  * (N_e / N_b)^2" , legend = false)

###
plot(les.t[20:end] ./ 86400, Δh[20:end] ./ h[20:end] .* (Nᵉ[20:end].^2) ./ N², ylims = (0.0, 0.5))
plot!(les.t[20:end] ./ 86400, Δh[20:end] ./ h1[20:end] .* 3.5, ylims = (0.05, 1.0))

#=
p2 = plot(les.t[20:end] ./ 86400, Δh[20:end] ./ h1[20:end], ylims = (0.0, 0.4)  )

plot(les.t[20:end] ./ 86400, Nᵉ[20:end].^2)
plot(les.t[20:end] ./ 86400, Ri[20:end])
=#
function moving_average(field, window)
    Nt = length(field)
    Nt_new = Nt - window
    low_pass_field = zeros(Nt_new)
    for i in 1:Nt_new
        low_pass_field[i] = mean(field[i:(i+window)])
    end
    return low_pass_field
end

###
Nt = length(les.t)
f_range = 50:Nt
ratio = Δh ./ h .* (Nᵉ.^2) ./ N²
dh = Δh ./ h
Cᴿ = Δh .* ΔB ./ (TKE_base )
nCᴿ = @. Nᵉ^2 * ( Δh^2 / ΔTKE)
# nCᴿ = @. Nᵉ^2 * ( Δh^2 / TKE_base)
low_pass_Ri = moving_average(Ri[f_range], 60)
low_pass_t = moving_average(les.t[f_range], 60)
low_pass_ratio = moving_average(ratio[f_range], 60)
low_pass_dh = moving_average(dh[f_range], 60)
low_pass_Cᴿ = moving_average(Cᴿ[f_range], 60)
low_pass_nCᴿ = moving_average(nCᴿ[f_range], 60)
low_pass_Nᵉ = moving_average(Nᵉ[f_range], 60)
###

plot(low_pass_t ./ 86400, 1 ./ low_pass_Ri .* 30)
plot(low_pass_t ./ 86400, low_pass_ratio, ylims = (0.2, 0.4), label = "dh/h x (N_e / N)^2", linewidth = 2)
# scatter(low_pass_Ri, low_pass_ratio, ylims = (0.2, 0.4))

plot!(low_pass_t ./ 86400, low_pass_dh .* 3, ylims = (0.2, 0.4), xlabel = "Days", label = "dh / h", ylabel = "Dimensionless Ratio" , legend = :topright, linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box)


###
histogram(low_pass_ratio, bins = 50, xlims = (0.2, 0.4))
# plot(les.t, dh, ylims = (0,0.2))

###

plot(les.t[20:end] ./ 86400, Cᴿ[20:end], ylims = (10^2,2 * 10^3), ylabel = "Richardson Number" )



plot(low_pass_t ./ 86400, low_pass_nCᴿ, ylims = (10^1,3 * 10^3), ylabel = "Low Pass Richardson Number" )


plot(low_pass_t ./ 86400, low_pass_Cᴿ, ylabel = "Low Pass Richardson Number", linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlabel = "days", legend = false , ylims = (1, 2))

plot(low_pass_t ./ 86400, low_pass_dh, ylabel = "dh/h", linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlabel = "days", legend = false , ylims = (0, 0.2), label = "dh/h")

plot(low_pass_t ./ 86400, low_pass_ratio, ylabel = "ratio", linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlabel = "days", legend = false , ylims = (0, 0.5))

plot(low_pass_t ./ 86400, (low_pass_Nᵉ).^2 ./ N², ylabel = "stratification ratio", linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box, xlabel = "days", legend = false , ylims = (1.5, 5.0))


plot(les.t, ΔB ./ Δh ./ (Nᵉ.^2), ylims = (0,1))


plot(les.t,  ΔB ./ 𝒮 ./ h2, ylims = (0,0.1))

plot(les.t ./86400 ,  (Nᵉ .^2) ./ 𝒮 , ylims = (0,0.5))


###
plot(les.t, sqrt.( ww_top ./ ww_star ), ylims = (0, 10^(-1)))
# plot(les.t, ww_base, ylims = (0, 10^(-7)))


###
field = les.T
max_field = maximum(field)
min_field = minimum(field)
animation = true
end_ind = floor(Int, Nt/2)
i = end_ind- 20
ϕ = field[:,i]
time_string = @sprintf("%.1i", les.t[i] ./ 86400)
p1 = plot(ϕ, les.z,  label = "Temperature", linewidth = 2, grid = true, gridstyle = :dash, gridalpha = 0.25, framestyle = :box , ylabel = "Depth [meters]", xlabel = "Temperature [Celcius]", xlims = (19.4, 19.6), ylims = (-60, 0), legend = false)
# title = "Day = " * time_string,


plot(les.t ./ 86400, Δh, ylims = (0,8))

plot(les.t ./ 86400, Δh .* (Nᵉ .^2) ./ (h .* N²), ylims = (0,0.5) )

plot(les.t ./ 86400, ΔB ./ (h .* N²), ylims = (0,0.3) )

plot(ΔB ./ Δh ./ (Nᵉ .^2))


plot(les.t, h, label = "measure 1")
plot!(les.t, h1, label = "measure 2")
plot!(les.t, h2, label = "measure 3")
analytic = @. ( 2.88 * les.t * Qᵇ / N² )^(0.5)
plot!(les.t, analytic, label = "analytic", legend = :topleft)

scatter(les.wT[:,end], les.z, legend = false)

minimum(les.wT[:,end]) ./ maximum(les.wT[:,end])
