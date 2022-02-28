# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "RBCL1" (RBC + linex cap. adj costs, first order linearization around SS),
# where LL is evaluated using BSPF



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Load external functions needed for this script:

include(sMyPath * "Functions/fModelFunctions.jl")
include(sMyPath * "Functions/fBSPF_linearSS.jl")


# Parameter labels (under this model) :

vParLabs    = ["rho_z","sig_z","rho_b","sig_b","r","tau","alpha","delta","nu","phi1","phi2","sigma_y","sigma_i","sigma_l"]


# Final function that computes LL under this model:

function fLL_L1bspf(θ)

      Mbspf = 400

      Ψ0, Ψ1, Ψu, Φ1, Φε, vYlabs, vXlabs = fGetStateSpaceRep(θ)

      return fBSPF_linearSS(mData,Mbspf,Ψ0,Ψ1,Φ1,Φε,1,Ψu*Ψu')[1]

end


# Only used for assessing BSPF accuracy (in "main_RBC_CheckBSPF.jl"):
function fLL_L1bspf(θ,Mbspf)

      Ψ0, Ψ1, Ψu, Φ1, Φε, vYlabs, vXlabs = fGetStateSpaceRep(θ)

      return fBSPF_linearSS(mData,Mbspf,Ψ0,Ψ1,Φ1,Φε,1,Ψu*Ψu')[1]

end
