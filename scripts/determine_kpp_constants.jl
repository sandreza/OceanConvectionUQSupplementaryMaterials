# utilizing notation of (Large 1994)
# OCEANIC VERTICAL MIXING: A REVIEW AND A MODEL WITH A NONLOCAL BOUNDARY LAYER PARAMETERIZATION
Cᵛ = 1.7    # should be between 1 and 2.1
βᵀ = -0.2   # entrainment law, page 367
Riᶜ = 0.3   # critical richardson number

cˢ = 98.96  # similarity constant, page 392
κ = 0.4     # von karman constant
ϵ = 0.1     # surface layer fraction, page 371

Cᴷᴱ = Cᵛ * (-βᵀ / (cˢ * ϵ))^(1/2) / (Riᶜ * κ) * (cˢ * ϵ * κ)^(1/3) #formula on page 372
# Cᴷᴱ = Cᵛ * (-βᵀ)^(1/2) / (Riᶜ * κ^(2/3)) * (cˢ * ϵ )^(-1/6)
# note that the critical richardson number drops out of the resulting expression in the strongly convective limit

# nonlocal diffusivity amplitude
Cstar = 10.0 # page 371
Cˢ = Cstar * κ * (cˢ * κ * ϵ)^(1/3) # page 371


# diffusivity amplitude
Cᴰ = κ * (cˢ * κ * ϵ )^(1/3)  #page 371

# in the paper we are using as the default
default_𝑪 = randn(4)

default_𝑪[1] = ϵ # the surface layer fraction
default_𝑪[2] = Cˢ
default_𝑪[3] = 1.36 # taken from ocean turb
default_𝑪[4] = Cᴷᴱ * Riᶜ  # a product of a bunch of constants. But only the Cke parameter enters in oceanturb
