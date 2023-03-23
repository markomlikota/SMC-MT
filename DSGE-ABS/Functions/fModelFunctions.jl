# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Various functions pertaining to the model (analytical model; RBC with cap. adj. costs) that can be used by several "models" (solution methods)



# -------------------------------------------------------------------------------


# Functions to go from θ to Model (SS-rep.) objects:

function fGetSteadyState(θ)

    # Extract parameters from θ:

    rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, ηdy, ηdwnom, ηr, ηinfl = θ

    β     = 1/(1+rA/400)
    πStar = 1 + πA/400
    γ     = 1 + γA/400
    φp    = (1-λpSS)/(λpSS*κ*πStar^2)


    gSS         = log(gStar)
    aSS         = log(γ)
    lambda_pSS  = log(λpSS)
    eSS         = 0
    wSS         = log( 1-λpSS )
    RSS         = log( γ/β * πStar )
    cSS         = log( (((1-λpSS)*(1-λwSS) * gStar^(-1/ν) )/χh )^(1/(τ+1/ν)) )
    ySS         = gSS + cSS
    dySS        = 0
    inflSS      = log(πStar)
    dwSS        = aSS
    dwnomSS     = dwSS + inflSS
    Phi_pSS     = 0
    Phi_wSS     = 0
    Phi_p_derSS = 0
    Phi_w_derSS = 0
    muSS        = λpSS


    return gSS, aSS, lambda_pSS, eSS, wSS, RSS, ySS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS

end
