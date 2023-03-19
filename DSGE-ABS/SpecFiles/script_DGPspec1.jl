# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# DGP specification file



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


### SAMPLE SIZE

T     = 92


### DEFINE TRUE VALUES FOR MODEL PARAMETERS


rA_0    = 1.88 # defines β     = 1/(1+rA/400)
πA_0    = 3.34 # defines πStar = 1 + πA/400
γA_0    = 1.98 # defines γ     = 1 + γA/400

τ_0     = 4.10
ν_0     = 0.10

κ_0     = 0.21 # defines φp    = (1-λpSS)/(λpSS*κ*πStar^2)
φw_0    = 11.7
ψw_0    = 59.4
ψp_0    = 165.0

ψ1_0    = 2.57
ψ2_0    = 0.79

ρr_0    = 0.73
ρg_0    = 0.96
ρa_0    = 0.07
ρp_0    = 0.9

σr_0    = 0.0017
σg_0    = 0.0083
σa_0    = 0.0047
σp_0    = 0.0654

gStar_0 = 1/0.85
λpSS_0  = 0.1
λwSS_0  = 0.1
χh_0    = 1


# StDs of measurement errors:
# Based on mL1:
# σ_dy_0 = 0.2013203131779747
# σ_dwnom_0 = 0.1565384665803447
# σ_r_0 = 0.2655211403716596
# σ_infl_0 = 0.1958441301967525

# Based on "true model" (VFI/L2/L3):
#L2:
σ_dy_0 = 0.2707567016393764
σ_dwnom_0 = 0.2073463519818686
σ_r_0 = 0.4203694560630014
σ_infl_0 = 0.3340350287056024
# #L3:
# σ_dy_0 = ...

# (values taken from simulating data based on corresponding DGP without MEs)



### DEFINE VECTOR OF PARAMETERS

θ0          = [rA_0, πA_0, γA_0, τ_0, ν_0, κ_0, φw_0, ψw_0, ψp_0, ψ1_0, ψ2_0, ρr_0, ρg_0, ρa_0, ρp_0, σr_0, σg_0, σa_0, σp_0, gStar_0, λpSS_0, λwSS_0, χh_0, σ_dy_0, σ_dwnom_0, σ_r_0, σ_infl_0]

nP          = length(θ0)

vParLabs    = ["rA","piA","gammaA","tau","nu","kappa","phi_w","psi_w","psi_p","psi1","psi2","rho_r","rho_g","rho_a","rho_p","sig_r","sig_g","sig_a","sig_p","gStar","lambda_p_ss","lambda_w_ss","chi_h","sig_eta_dy","sig_eta_dwnom","sig_eta_r","sig_eta_infl"]
