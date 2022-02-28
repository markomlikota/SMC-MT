# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# DGP specification file



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


### SAMPLE SIZE

T     = 80


### DEFINE TRUE VALUES FOR MODEL PARAMETERS


# Parameters for exogenous processes:

ρz_0       = 0.95
σz_0       = 0.0000001 #0.005

ρb_0       = 0.90
σb_0       = 0.0000001 #0.004


# SS parameters:

r_0        = 3
τ_0        = 2.0
α_0        = 0.35
δ_0        = 0.08
ν_0        = 1.0


# # Back out Zstar and Bstar s.t. ySS = lSS = 1:
#
# β_0        = 1/(1 + r_0/400)
#
# Zstar_0    = ( (1/β_0 - (1-δ_0))/α_0 )^α_0
# Bstar_0    = (1-α_0) * ( 1 - δ_0*α_0/(1/β_0 - (1-δ_0)) )^(-τ_0)


# Capital adjustment cost parameters:

ϕ1_0       = 0
ϕ2_0       = 0


# StDs of measurement errors:
# Based on mL1:
# σy_0 = 0.0032926223722702
# σi_0 = 0.0056619204952717
# σl_0 = 0.0020552117164827
# Based on mVFI:
σy_0 = 0.0035464492908291
σi_0 = 0.0064019169274366
σl_0 = 0.0010346776618169
# (values taken from simulating data based on corresponding DGP without MEs)


### DEFINE VECTOR OF PARAMETERS

θ0          = [ρz_0, σz_0, ρb_0, σb_0, r_0, τ_0, α_0, δ_0, ν_0, ϕ1_0, ϕ2_0, σy_0, σi_0, σl_0]

nP          = length(θ0)

vParLabs    = ["rho_z","sig_z","rho_b","sig_b","r","tau","alpha","delta","nu","phi1","phi2","sigma_y","sigma_i","sigma_l"]
