# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# SMC specification file:



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# SMC tuning parameters:

N               = 1000
useAdaptive     = 1
λ               = 2
Nϕ              = 200
# ϕ_Nϕ            = 1.00 # value for phi in last stage (fixed TS can target it exactly; adaptive TS terminates when phi becomes larger or equal to ϕ_Nϕ)
α               = 0.95
Nmh             = 2
nB              = 1 # number of blocks for RWMH in mutation steps, random blocking is used
c0              = 0.5
accStar         = 0.25
Nbar            = N/2
showMessages    = 0


# Note: once nP (length of θ) and ϕ_Nϕ are defined, the above will be combined into tSettings, which is argument of fSMC():
# tSettings       = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages) #(use tuple to keep different types (Int vs Float))
