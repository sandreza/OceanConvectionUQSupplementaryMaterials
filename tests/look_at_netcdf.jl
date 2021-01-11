using NCDatasets, Statistics, Plots
files = [
    pwd() * "/three_layer_constant_fluxes_cubic_hr48_Qu0.0e+00_Qb1.0e-08_f1.0e-04_Nh256_Nz128_free_convection_2days_Qb1e-8.nc",
    pwd() * "/checkthis.nc"
]
file = files[2]
ds = Dataset(file,"r")
field = ds["T"]
t = ds["time"][:]
meanT = zeros(length(t))
key = keys(field)
for i in eachindex(t) 
    meanT[i] = mean(field[:,i])  
end

plot(t, meanT)

αg = 0.00196133
for  i in 2:length(t)
    println(αg * (meanT[i-1] - meanT[i]) / (t[i] - t[i-1]))
end