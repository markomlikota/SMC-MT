# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# May 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "L1" (first order linearization around SS)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Load external functions needed for this script:

include(sMyPath * "Functions/fModelFunctions.jl")
include(sMyPath * "Functions/fKFLL.jl")


# Parameter labels (under this model) :

vParLabs    = ["rA","piA","gammaA","tau","nu","kappa","phi_w","psi_w","psi_p","psi1","psi2","rho_r","rho_g","rho_a","rho_p","sig_r","sig_g","sig_a","sig_p","gStar","lambda_p_ss","lambda_w_ss","chi_h","sig_eta_dy","sig_eta_dwnom","sig_eta_r","sig_eta_infl"]


# Model Solution and Particle Filter at given θ:

myInv(x) = inv(x) #pinv(x; atol=1e-15)

#Load pre-processed model:
sPathModelspecSolveDSGE_processed      = sMyPath * "SpecFiles/" * "script_ModelspecLogLX_SolveDSGE_processed.txt"
modelStructureSolveDSGE                = retrieve_processed_model(sPathModelspecSolveDSGE_processed)



function fSolveModelL1(θ)


      # Assign parameters to model pre-processed by SolveDSGE:

      rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, σ_dy, σ_dwnom, σ_r, σ_infl = θ

      ψw == 0 ? ψw = 1e-5 : nothing
      ψp == 0 ? ψp = 1e-5 : nothing
      β     = 1/(1+rA/400)
      πStar = 1 + πA/400
      γ     = 1 + γA/400
      φp    = (1-λpSS)/(λpSS*κ*πStar^2)
      vParams_SolveDSGE             = [τ, β, λpSS, λwSS, φw, ψw, φp, ψp, ρa, σa, γ, ρg, σg, gStar, ρr, ψ1, ψ2, σr, πStar, χh, ν, ρp, σp]
      modelStructureAtThisθ       = assign_parameters(modelStructureSolveDSGE,vParams_SolveDSGE)


      # Compute 1st order linearization around SS:

      gSS, aSS, lambda_pSS, eSS, wSS, RSS, ySS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS = fGetSteadyState(θ)

      vSS_SolveDSGE                 = [gSS, aSS, lambda_pSS, eSS, ySS, RSS, wSS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS] #has to be in order specified by SolveDSGE in modelStructureAtThisθ

      modelL1_SolveDSGE             = PerturbationScheme(vSS_SolveDSGE,1.0,"first")
      modelL1_SolveDSGE_sol         = solve_model(modelStructureAtThisθ,modelL1_SolveDSGE)


      # Put in linear SS form:

      Φ1    = modelL1_SolveDSGE_sol.hx
      Φε    = modelL1_SolveDSGE_sol.k

      Ψ0    = [γA/4, γA/4+πA/4, log(γ*πStar/β)*400,πA]
      Ψ1    = zeros(4,7)
      Ψ1[1,:] = [0,100,zeros(8)...]'*modelL1_SolveDSGE_sol.gx + [0,100,zeros(5)...]' # 100 * (dy + a)
      Ψ1[2,:] = [0,0,100,zeros(7)...]'*modelL1_SolveDSGE_sol.gx                      # 100 * dwnom
      Ψ1[3,:] = [zeros(5)...,400,0]                                                  # 400 * R
      Ψ1[4,:] = [0,0,0,400,zeros(6)...]'*modelL1_SolveDSGE_sol.gx                    # 400 * infl

      Σ_u   = diagm([σ_dy, σ_dwnom, σ_r, σ_infl].^2)

      return Ψ0,Ψ1,Φ1,Φε,Σ_u

end

# Log-likelihood:

function fLL_L1(θ)

      Ψ0, Ψ1, Φ1, Φε, Σ_u = fSolveModelL1(θ)

      vLL = fKalmanFilter(mData,Ψ0,Ψ1,Φ1,Φε,0,1,Σ_u)[4]

      return sum(vLL)

end



function fSimulate_mL1(θ,T,randomSeedNumber,randomSeedNumber2)


      # Assign parameters to model pre-processed by SolveDSGE:

      rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, σ_dy, σ_dwnom, σ_r, σ_infl = θ

      ψw == 0 ? ψw = 1e-5 : nothing
      ψp == 0 ? ψp = 1e-5 : nothing
      β     = 1/(1+rA/400)
      πStar = 1 + πA/400
      γ     = 1 + γA/400
      φp    = (1-λpSS)/(λpSS*κ*πStar^2)
      vParams_SolveDSGE             = [τ, β, λpSS, λwSS, φw, ψw, φp, ψp, ρa, σa, γ, ρg, σg, gStar, ρr, ψ1, ψ2, σr, πStar, χh, ν, ρp, σp]
      modelStructureAtThisθ       = assign_parameters(modelStructureSolveDSGE,vParams_SolveDSGE)


      # Compute 1st order linearization around SS:

      gSS, aSS, lambda_pSS, eSS, wSS, RSS, ySS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS = fGetSteadyState(θ)

      vSS_SolveDSGE                 = [gSS, aSS, lambda_pSS, eSS, ySS, RSS, wSS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS] #has to be in order specified by SolveDSGE in modelStructureAtThisθ

      modelL1_SolveDSGE             = PerturbationScheme(vSS_SolveDSGE,1.0,"first")
      modelL1_SolveDSGE_sol         = solve_model(modelStructureAtThisθ,modelL1_SolveDSGE)


      mAllVars        = SolveDSGE.simulate(modelL1_SolveDSGE_sol,[gSS, aSS, lambda_pSS, eSS, ySS, RSS, wSS],T; rndseed=randomSeedNumber)

      mX      = mAllVars[1:7,:]'
      mY      = mAllVars[8:end,:]'

      mObservables = zeros(T,4)
      mObservables[:,1] = γA/4 .+ 100*(mX[:,2].-dySS + mY[:,2].-aSS)
      mObservables[:,2] = γA/4+πA/4 .+ 100*(mY[:,3].-dwnomSS)
      mObservables[:,3] = log(γ*πStar/β)*400 .+ 400*(mX[:,6].-RSS)
      mObservables[:,4] = πA .+ 400*(mY[:,4].-inflSS)
      #modelStructureSolveDSGE.variables


      Random.seed!(randomSeedNumber2)

      mu = diagm([σ_dy, σ_dwnom, σ_r, σ_infl]) * randn(4,T)
      mObservables .+= mu'

      return mObservables

end
