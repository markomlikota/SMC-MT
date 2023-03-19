# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# May 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "L2" (second-order lin. around SS)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# using SolveDSGE



# Load external functions needed for this script:

include(sMyPath * "Functions/fModelFunctions.jl")
include(sMyPath * "Functions/fMyMNDraw.jl")


# Parameter labels (under this model) :

vParLabs    = ["rA","piA","gammaA","tau","nu","kappa","phi_w","psi_w","psi_p","psi1","psi2","rho_r","rho_g","rho_a","rho_p","sig_r","sig_g","sig_a","sig_p","gStar","lambda_p_ss","lambda_w_ss","chi_h","sig_eta_dy","sig_eta_dwnom","sig_eta_r","sig_eta_infl"]


# Model Solution and Particle Filter at given θ:

myInv(x) = inv(x) #pinv(x; atol=1e-15)

#Load pre-processed model:
sPathModelspecSolveDSGE_processed      = sMyPath * "SpecFiles/" * "script_ModelspecLogLX_SolveDSGE_processed.txt"
modelStructureSolveDSGE                = retrieve_processed_model(sPathModelspecSolveDSGE_processed)



function fSolveModelL2(θ)


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


      # Compute 2nd order linearization around SS:

      gSS, aSS, lambda_pSS, eSS, wSS, RSS, ySS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS = fGetSteadyState(θ)

      vSS_SolveDSGE                 = [gSS, aSS, lambda_pSS, eSS, ySS, RSS, wSS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS] #has to be in order specified by SolveDSGE in modelStructureAtThisθ

      modelL2_SolveDSGE             = PerturbationScheme(vSS_SolveDSGE,1.0,"second")
      modelL2_SolveDSGE_sol         = solve_model(modelStructureAtThisθ,modelL2_SolveDSGE)


      # Obtain function to simulate forward _vs to vs and obtain the corresponding vector of observables:

      Ψ0    = [γA/4, γA/4+πA/4, log(γ*πStar/β)*400,πA]
      function fΦΨ(_vs)
          mVars = SolveDSGE.simulate(modelL2_SolveDSGE_sol,_vs,2;seed=Int(floor(rand()*1e10)))
          vs = mVars[1:7,2]
          vy = mVars[8:end,2]

          a = vs[2]
          R = vs[6]
          dy = vy[2]
          dwnom = vy[3]
          infl = vy[4]

          vObservables = Ψ0 + [100*(dy-dySS + a-aSS), 100*(dwnom-dwnomSS), 400*(R-RSS), 400*(infl-inflSS)]

          return vs, vObservables
      end


      # Obtain function to initialize BSPF based on 1st order lin.:

      modelL2_SolveDSGE_L1             = PerturbationScheme(vSS_SolveDSGE,1.0,"first")
      modelL2_SolveDSGE_L1_sol         = solve_model(modelStructureAtThisθ,modelL2_SolveDSGE_L1)

      Φ1_SolveDSGE  = modelL2_SolveDSGE_L1_sol.hx
      Φε_SolveDSGE  = modelL2_SolveDSGE_L1_sol.k
      nSlSS         = size(Φ1_SolveDSGE,1)

      P00vec      = myInv( diagm(vec(ones(nSlSS^2,1))) - kron(Φ1_SolveDSGE,Φ1_SolveDSGE) ) * vec(Φε_SolveDSGE * Φε_SolveDSGE')
      P00         = reshape(P00vec,(nSlSS,nSlSS))

      function fDrawInitialStates()
          vs0_SSdev       = fMyMNDraw(zeros(nSlSS),P00)
          vs0             = [gSS, aSS, lambda_pSS, eSS, ySS, RSS, wSS] .+ vs0_SSdev
          return vs0
      end

      return fΦΨ, fDrawInitialStates

end


ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)
# #apparently, this should be faster: (BUT: CHECK WHETHER IT GIVES SAME RESULT AS ABOVE)
# ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * ( Σ\(x-μ) )


function fBSPF(mData,Mbspf,nS,θ,fΦΨ,fDrawInitialStates)

    # Bootstrap Particle Filter for linear state space of the form
    #           s_t  =  fΦ(s_{t-1},ε_t;θ)         ,   ε_t ∼ fFε(ε_{t-1};θ)
    #           y_t  =  fΨ(s_t;θ)          +  u_t ,   u_t ∼ N(0,Σ_u(θ))
    #
    # Inputs:
    # - (T x nY) matrix of data for likelihood evaluation,
    # - number of particles, number of endogenous states
    # - vector of model parameters
    # - SS function (which encodes the model solution at this particular θ and simulates s_{t-1} forward to obtain s_t and y_t)
    # - function to draw initial states (s_0) (which also was generated based on the model solution at this particular θ)
    #
    # Output: loglik, array of particles, matrix of weights
    #


    ### ------------------------ SET UP ----------------------------------------

    T,nY        = size(mData)

    rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, σ_dy, σ_dwnom, σ_r, σ_infl = θ

    Σ_u         = diagm([σ_dy, σ_dwnom, σ_r, σ_infl].^2)


    ### ------------------------ INITIALIZE OBJECTS ----------------------------

    aParticlesEndoStates    = zeros(T+1,Mbspf,nS)
    mWeights                = zeros(T+1,Mbspf)
    vLogLik                 = zeros(T)



    ### ------------------------ INITIALIZE ALGORITHM (RECURSIONS) -------------

    [aParticlesEndoStates[1,jj,:] = fDrawInitialStates() for jj = 1:Mbspf]

    mWeights[1,:]       .= 1


    ### ------------------------ RECURSIONS ------------------------------------

    for t = 1:T

        # --- A: Forecasting s_t and y_t : -------------------------------------------------
        # compute s_t = fΦ(s_{t-1},ε_t;θ) and y_t = fΨ(s_t;θ)

        vv_s_Obs                     = [ fΦΨ(aParticlesEndoStates[t,jj,:]) for jj = 1:Mbspf ]


        # evaluate p(y_t|s_t,θ) (coming from y_t = fΨ(s_t;θ) + u_t ,   u_t ∼ N(0,Σ_u(θ)) )

        vLogIncrWeights         = [ ll( mData[t,:], vv_s_Obs[jj][2] , Σ_u )  for jj = 1:Mbspf ]


        # Likelihood increment:

        vHelp                   = vLogIncrWeights .+ log.(mWeights[t,:]) #log of w-tilde * W
        c                       = maximum(vHelp)
        vHelp                   = vHelp .- c

        vLogLik[t]              = log(mean(exp.(vHelp))) + c


        # --- B: Normalized Weights : ----------------------------------------------

        vNormWeights            = exp.(vHelp) ./ mean(exp.(vHelp))


        # --- C: Selection : --------------------------------------------------------

        ESS                     = Mbspf / mean( vNormWeights.^2 )
        doSelect                = ESS < Mbspf/2

        if doSelect == 1

            vSelectionIndices       = resampleNYFED(vNormWeights)
            [aParticlesEndoStates[t+1,jj,:]   = vv_s_Obs[vSelectionIndices[jj]][1] for jj = 1:Mbspf]
            mWeights[t+1,:]         .= 1

        elseif doSelect == 0

            mWeights[t+1,:]         = vNormWeights
            [aParticlesEndoStates[t+1,jj,:]   = vv_s_Obs[jj][1] for jj = 1:Mbspf]

        end

    end


    ### ------------------------ COMPUTE LL APPROX. ----------------------------

    loglik          = sum( vLogLik )


    return loglik, aParticlesEndoStates, mWeights, vLogLik

end



# -- Resulting Likelihood Function ---------------------------------------------


function fLL_L2(θ)

    fΦΨ, fDrawInitialStates = fSolveModelL2(θ)

    Mbspf   = 25000

    nS      = 7
    nε      = 4

    return fBSPF(mData,Mbspf,nS,θ,fΦΨ,fDrawInitialStates)[1]

end


# Only used for assessing BSPF accuracy (in "main_RBC_CheckBSPF.jl"):


function fLL_L2(θ,Mbspf)

    timeSolution = @elapsed fΦΨ, fDrawInitialStates = fSolveModelL2(θ)

    nS      = 7
    nε      = 4

    timePF = @elapsed loglik = fBSPF(mData,Mbspf,nS,θ,fΦΨ,fDrawInitialStates)[1]

    return loglik, timeSolution, timePF

end



# -- Function to simulate data based on this nonlinear SS rep: -----------------

function fSimulate_mL2(θ,T,randomSeedNumber,randomSeedNumber2)

    rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, σ_dy, σ_dwnom, σ_r, σ_infl = θ

    ψw == 0 ? ψw = 1e-5 : nothing
    ψp == 0 ? ψp = 1e-5 : nothing
    β     = 1/(1+rA/400)
    πStar = 1 + πA/400
    γ     = 1 + γA/400
    φp    = (1-λpSS)/(λpSS*κ*πStar^2)
    vParams_SolveDSGE             = [τ, β, λpSS, λwSS, φw, ψw, φp, ψp, ρa, σa, γ, ρg, σg, gStar, ρr, ψ1, ψ2, σr, πStar, χh, ν, ρp, σp]
    modelStructureAtThisθ       = assign_parameters(modelStructureSolveDSGE,vParams_SolveDSGE)


    # Compute 2nd order linearization around SS:

    gSS, aSS, lambda_pSS, eSS, wSS, RSS, ySS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS = fGetSteadyState(θ)

    vSS_SolveDSGE                 = [gSS, aSS, lambda_pSS, eSS, ySS, RSS, wSS, cSS, dySS, dwnomSS, inflSS, dwSS, Phi_pSS, Phi_wSS, Phi_p_derSS, Phi_w_derSS, muSS] #has to be in order specified by SolveDSGE in modelStructureAtThisθ

    modelL2_SolveDSGE             = PerturbationScheme(vSS_SolveDSGE,1.0,"second")
    modelL2_SolveDSGE_sol         = solve_model(modelStructureAtThisθ,modelL2_SolveDSGE)


    # Simulate data:

    mAllVars        = SolveDSGE.simulate(modelL2_SolveDSGE_sol,[gSS, aSS, lambda_pSS, eSS, ySS, RSS, wSS],T; rndseed=randomSeedNumber)

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
