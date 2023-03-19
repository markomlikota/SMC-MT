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

# # with variables not in logs:
# function fGetSteadyState(θ)
#
#     # Extract parameters from θ:
#
#     rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, ηdy, ηdwnom, ηr, ηinfl = θ
#
#     β     = 1/(1+rA/400)
#     πStar = 1 + πA/400
#     γ     = 1 + γA/400
#     φp    = (1-λpSS)/(λpSS*κ*πStar^2)
#
#
#     # Fix ySS, lSS:
#
#     gSS         = log(gStar)
#     aSS         = log(γ)
#     lambda_pSS  = log(λpSS)
#     eSS         = 0
#     wSS         = log( 1-λpSS )
#     RSS         = γ/β * πStar
#     cSS         = (((1-λpSS)*(1-λwSS) * gStar^(-1/ν) )/χh )^(1/(τ+1/ν))
#     ySS         = gSS + cSS
#     dySS        = 0
#     inflSS      = log(πStar)
#     dwSS        = aSS
#     dwnomSS     = dwSS + inflSS
#     Phi_pSS     = 0
#     Phi_wSS     = 0
#     Phi_p_derSS = 0
#     Phi_w_derSS = 0
#     muSS        = λpSS
#
#
#     return gSS, aSS, lambda_pSS, eSS, wSS, RSS, ySS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS
#
# end


# function fGetStateSpaceRep(θ)
#
#     # Given a parameter vector, this function computes the system matrices for the state space representation. They are of the following form:
#     #   Measurement Equation:   y_t = Ψ0 + Ψ1 s_t     + Ψu u_t
#     #   Transition Equation:    s_t =      Φ1 s_{t-1} + Φε ε_t
#
#     # Thereby, it uses the procedure in Sims(2002) to solve the rational expectation system of equations (Schur-decomposition)
#
#
#     ### ------------------------ EXTRACT PARAMETERS FROM Θ ---------------------
#
#     ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ
#
#
#     # Define "helping" parameters and others:
#
#     β       = 1/(1 + r/400)
#
#     h1      = 1-β*(1-δ)
#     h2      = (1-α)/(α*ν+1)
#     h3      = ν/(α*ν+1)
#
#     cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)
#
#
#     ### -------- CREATE MATRICES FOR LINEAR RATIONAL EXPECTATIONS FORM ---------
#     #                 (  Γ0 s_t = Γ1 s_{t-1} + Ψ ε_t + Π η_t  )
#
#     col_y      = 1
#     col_c      = 2
#     col_i      = 3
#     col_l      = 4
#     col_k      = 5
#     col_z      = 6
#     col_b      = 7
#     col_Ec     = 8
#     col_Ek     = 9
#
#     vXlabs     = ["y","c","i","l","k'","z","b","Ec","Ek''"]
#
#     nX         = 9
#     mX         = 2 # how many forward-looking variables (for solving R.E. system below)
#
#     Γ0         = zeros(nX,nX)
#     Γ1         = zeros(nX,nX)
#     Ψ          = zeros(nX,2)
#     Π          = zeros(nX,2)
#
#     Σε         = [σz 0; 0 σb]
#
#
#     # Fill-in Γ0 & Γ1:
#
#     # equation for output (production function)
#
#     Γ0[col_y,col_y] = 1
#     Γ0[col_y,col_l] = -(1-α)
#     Γ0[col_y,col_z] = -1
#     Γ1[col_y,col_k] = α
#
#     # equation for consumption (goods market clearing)
#
#     Γ0[col_c,[col_y,col_c,col_i]] = [ -1, cSS/ySS, iSS/ySS ]
#
#     # equation for investment
#
#     Γ0[col_i,[col_i,col_k]] = [ δ, -1 ]
#     Γ1[col_i,col_k]         = -(1-δ)
#
#     # equation for labor supply:
#
#     Γ0[col_l,[col_c,col_l,col_b,col_y]] = [ τ, 1+1/ν, 1, -1 ]
#
#     # equation for savings/investment/capital tmrw choice (EE):
#
#     Γ0[col_k,col_c]    = -τ
#     Γ0[col_k,col_k]    = ϕ1*(1+β) + h1*h2
#     Γ0[col_k,col_b]    = h1*h2*ν*ρb
#     Γ0[col_k,col_z]    = -h1*(1+h2*ν)*ρz
#     Γ0[col_k,col_Ec]   = τ*(1+ν*h1*h2)
#     Γ0[col_k,col_Ek]   = -β*ϕ1
#
#     Γ1[col_k,col_k]    = ϕ1
#
#     # equation for productivity process z:
#
#     Γ0[col_z,col_z]    = 1
#     Γ1[col_z,col_z]    = ρz
#     Ψ[col_z,1]         = 1
#
#     # equation for taste process b:
#
#     Γ0[col_b,col_b]    = 1
#     Γ1[col_b,col_b]    = ρb
#     Ψ[col_b,2]         = 1
#
#     # equation for E[c']:
#
#     Γ0[col_Ec,col_c]   = 1
#     Γ1[col_Ec,col_Ec]  = 1
#     Π[col_Ec,1]        = 1
#
#     # equation for E[k'']:
#
#     Γ0[col_Ek,col_k]   = 1
#     Γ1[col_Ek,col_Ek]  = 1
#     Π[col_Ek,2]        = 1
#
#
#
#     ### -------- SOLVE LINEAR RATIONAL EXPECTATIONS SYSTEM ---------------------
#     # to obtain transition equation, using Schur decomposition (Schorfheide & Herbst (2015, p. 16ff.), Sims (2002)):
#
#     # Schur decomposition:
#     F        = schur(Γ0, Γ1)
#
#     Λ        = F.S
#     Ω        = F.T
#     Q        = F.Q'
#     Z        = F.Z
#     # isapprox(Q' * Λ * Z',Γ0)
#     # isapprox(Q' * Ω * Z',Γ1)
#
#
#     # Sort according so that largest two generalized Eigenvalues appear last:
#     vGE      = diag(Ω)./diag(Λ) #vector of generalized Eigenvalues
#
#     # select = abs.(vGE) .< 1
#     vIndices = sortperm(abs.(vGE))
#     vSelect  = vec(repeat([true],nX,1))
#     vSelect[vIndices[end-1:end]] .= false
#
#     ordschur!(F,vSelect)
#
#
#     v1        = 1:(nX-mX)
#     v2        = (nX-mX+1):nX
#     Z1  = Z[v1,:];    Z2  = Z[v2,:]
#     Q1  = Q[v1,:];    Q2  = Q[v2,:]
#     Λ11 = Λ[v1,v1];  Λ12 = Λ[v1,v2]; Λ21 = Λ[v2,v1]; Λ22 = Λ[v2,v2]
#     Ω11 = Ω[v1,v1];  Ω12 = Ω[v1,v2]; Ω21 = Ω[v2,v1]; Ω22 = Ω[v2,v2]
#
#     if mX == 1
#         Φ       = ( Q1* Π ) /( (Q2* Π )[1])
#     else
#         Φ       = ( Q1* Π ) * inv( Q2* Π )
#     end
#
#     m1      = [Λ11 Λ12-Φ*Λ22; zeros(mX,nX-mX) diagm(vec(ones(mX,1)))]*Z'
#     m2      = [Ω11 Ω12-Φ*Ω22; zeros(mX,nX) ]*Z'
#     m3      = [Q1-Φ*Q2 ; zeros(mX,nX)] * Ψ
#
#     Φ1      = ( inv(m1) * m2 )
#     Φε      = ( inv(m1) * m3 )
#
#     Φε      *= Σε
#
#
#
#     ### -------- COMPUTE MATRICES FOR MEASUREMENT EQUATION ---------------------
#
#     Ψ0          = zeros(3,1)
#     Ψ0[2]       = log(iSS)
#
#     Ψ1          = zeros(3,nX)
#     Ψ1[1,col_y] = 1
#     Ψ1[2,col_i] = 1
#     Ψ1[3,col_l] = 1
#
#     Ψu          = diagm([ σy, σi, σl ])
#
#     vYlabs      = ["log y^o","log i^o","log l^o"]
#
#
#     return Ψ0, Ψ1, Ψu, Φ1, Φε, vYlabs, vXlabs
#
# end
#
#
# function fProcessModelL2(sMyPath)
#
#     # Function to process ModelL2 using package "SolveDSGE"
#
#     sFileNameModelspecL2    = sMyPath * "SpecFiles/" * "script_ModelspecL2.txt"
#     sPathModelspecL2        = joinpath(@__DIR__,sFileNameModelspecL2)
#
#     process_model(sPathModelspecL2)
#
#     sFileNameModelspecL2_processed      = sMyPath * "SpecFiles/" * "script_ModelspecL2_processed.txt"
#     sFPathModelspecL2_processed         =  joinpath(@__DIR__,sFileNameModelspecL2_processed)
#
#     modelStructureL2        = retrieve_processed_model(sFPathModelspecL2_processed)
#
#     return modelStructureL2
#
# end
