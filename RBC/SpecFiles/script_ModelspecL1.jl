# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "RBCL1" (RBC + linex cap. adj costs, first order linearization around SS)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Load external functions needed for this script:

include(sMyPath * "Functions/fModelFunctions.jl")
include(sMyPath * "Functions/fKFLL.jl")


# Parameter labels (under this model) :

vParLabs    = ["rho_z","sig_z","rho_b","sig_b","r","tau","alpha","delta","nu","phi1","phi2","sigma_y","sigma_i","sigma_l"]


# Log-likelihood:

function fLL_L1(θ)

      Ψ0, Ψ1, Ψu, Φ1, Φε, vYlabs, vXlabs = fGetStateSpaceRep(θ)

      vLL = fKalmanFilter(mData,Ψ0,Ψ1,Φ1,Φε,0,1,Ψu*Ψu')[4]

      return sum(vLL)

end
