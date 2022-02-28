# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Jan 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "L3" (RBC + linex cap. adj costs, first order linearization around SS)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# using SolveDSGE



# Load external functions needed for this script:

include(sMyPath * "Functions/fModelFunctions.jl")
include(sMyPath * "Functions/fMyMNDraw.jl")


# Parameter labels (under this model) :

vParLabs    = ["rho_z","sig_z","rho_b","sig_b","r","tau","alpha","delta","nu","phi1","phi2","sigma_y","sigma_i","sigma_l"]


# Model Solution and Particle Filter at given θ:

myInv(x) = inv(x) #pinv(x; atol=1e-15)

#Load pre-processed model:
sPathModelspecSolveDSGE_processed      = sMyPath * "SpecFiles/" * "script_ModelspecLX_SolveDSGE_processed.txt"
modelStructureSolveDSGE                = retrieve_processed_model(sPathModelspecSolveDSGE_processed)

function fSolveModelL3(θ)


      cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)


      # Assign parameters to model pre-processed by SolveDSGE:

      ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

      ϕ2 == 0 ? ϕ2 = 1e-5 : nothing
      vParams_SolveDSGE             = [ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, Zstar, Bstar]
      modelStructureAtThisθ         = assign_parameters(modelStructureSolveDSGE,vParams_SolveDSGE)


      # Compute 2nd order linearization around SS:

      vSS_SolveDSGE                 = [0.0, 0.0, kSS, ySS, cSS, lSS, iSS, 0.0, 0.0]

      modelL3_SolveDSGE             = PerturbationScheme(vSS_SolveDSGE,1.0,"third")
      modelL3_SolveDSGE_sol         = solve_model(modelStructureAtThisθ,modelL3_SolveDSGE)


      # Put in State Space (SS) form and obtain functions:

      modelL3_SolveDSGE_sol_SS      = state_space_eqm(modelL3_SolveDSGE_sol)

      f_dec_rule                      = modelL3_SolveDSGE_sol_SS.g
      f_state_trans                   = modelL3_SolveDSGE_sol_SS.h

      # vy = dec_rule([0.00,0.00,kSS]) # inputs are (here) [z,b,k] (the state variables, with the order given by dsge.variables), output is [y,c,l,i,ac,acdev]
      # vs_ = state_trans([0.00,0.00,kSS],[0.00, 0.00]) # inputs are again [z,b,k] and also [eps] (the vector of shocks), output is [z',b',k']

      function fΦ(_vs,vε)
            # vε are shocks (εz,εb), each ∼N(0,1)
            # vs_ = [z,b,k], vs = [z',b',k']
            vs = f_state_trans(_vs,vε)
            return vs
      end

      function fΨ(vs)
            vy = f_dec_rule(vs)
            y,c,l,i,ac,acdev = vy
            return log.( [y,i,l] )
      end


      # Compute 1st order linearization to get initial states for PF:

      modelL3_SolveDSGE_L1             = PerturbationScheme(vSS_SolveDSGE,1.0,"first")
      modelL3_SolveDSGE_L1_sol         = solve_model(modelStructureAtThisθ,modelL3_SolveDSGE_L1)

      Φ1_SolveDSGE = modelL3_SolveDSGE_L1_sol.hx
      Φε_SolveDSGE = modelL3_SolveDSGE_L1_sol.k
      nSlSS = size(Φ1_SolveDSGE,1)

      P00vec      = myInv( diagm(vec(ones(nSlSS^2,1))) - kron(Φ1_SolveDSGE,Φ1_SolveDSGE) ) * vec(Φε_SolveDSGE * Φε_SolveDSGE')
      P00         = reshape(P00vec,(nSlSS,nSlSS))

      function fDrawInitialStatesBSPF()
          vs0_SSdev       = fMyMNDraw(zeros(nSlSS),P00)
          z,b,khat        = vs0_SSdev
          k               = kSS*(1+khat)
          vs0             = [z,b,k]
          vε0             = randn(2)
          return vs0, vε0
      end


      return fΦ, fΨ, fDrawInitialStatesBSPF

end


ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)
# #apparently, this should be faster: (BUT: CHECK WHETHER IT GIVES SAME RESULT AS ABOVE)
# ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * ( Σ\(x-μ) )


function fBSPF(mData,Mbspf,nS,nε,θ,fΨ,fΦ,fDrawInitialStates)

    # Bootstrap Particle Filter for linear state space of the form
    #           s_t  =  fΦ(s_{t-1},ε_t;θ)         ,   ε_t ∼ fFε(ε_{t-1};θ)
    #           y_t  =  fΨ(s_t;θ)          +  u_t ,   u_t ∼ N(0,Σ_u(θ))
    # Note: with slight abuse of notation, ε are not the shocks/innovations to the exogenous shock processes, but the exogenous processes themselves (e.g. TFP process z that evolves acc. to AR(1))
    #
    # Inputs:
    # - (T x nY) matrix of data for likelihood evaluation,
    # - number of particles, number of endogenous states, number of exogenous states
    # - vector of model parameters
    # - SS objects (functions) (which encode the model solution at this particular θ)
    # - function to draw initial states (vs,vε) (which also was generated based on the model solution at this particular θ)
    #
    # Output: loglik, array of particles, matrix of weights
    #


    ### ------------------------ SET UP ----------------------------------------

    T,nY        = size(mData)

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    # cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)

    Σ_u         = diagm([σy,σi,σl].^2)


    ### ------------------------ INITIALIZE OBJECTS ----------------------------

    aParticlesEndoStates    = zeros(T+1,Mbspf,nS)
    aParticlesExoStates     = zeros(T+1,Mbspf,nε)
    mWeights                = zeros(T+1,Mbspf)
    vLogLik                 = zeros(T)



    ### ------------------------ INITIALIZE ALGORITHM (RECURSIONS) -------------

    for jj = 1:Mbspf

        vs, vε      = fDrawInitialStates()

        aParticlesEndoStates[1,jj,:] = vs

        aParticlesExoStates[1,jj,:] = vε

        # Note that in this particular case, because of the way the model solution is obtained (& SS-rep. functions are written),
        # ε_t is [εz, εb] rather than [z,b], so it wouldn't need to be drawn initially
        # Still, this is a tiny inefficiency, and I kept the code structure general

    end

    mWeights[1,:]       .= 1


    # vVarHere = aParticlesEndoStates[1,:,4]
    # histogram(vVarHere)

    # vVarHere = aParticlesExoStates[1,:,2]
    # histogram(vVarHere)



    ### ------------------------ RECURSIONS ------------------------------------

    for t = 1:T


        # --- a1) Forecasting ε_t : -------------------------------------------------
        # draw from p(ε_t|ε_{t-1},θ) = fFε(ε_{t-1},θ)

        vvε                 = [ randn(2) for jj = 1:Mbspf ]

        # Note that in this particular case, because of the way the model solution is obtained (& SS-rep. functions are written),
        # ε_t is [εz, εb] rather than [z,b], so it doesn't depend on ε_{t-1}

        # vVarHere = getindex.(vvε,2)
        # histogram(vVarHere)


        # --- a2) Forecasting s_t : -------------------------------------------------
        # compute s_t = fΦ(s_{t-1},ε_t;θ)

        vvs                     = [ fΦ(aParticlesEndoStates[t,jj,:],vvε[jj]) for jj = 1:Mbspf ]

        # vVarHere = getindex.(vvs,4)
        # histogram(vVarHere)


        # --- b) Forecasting y_t : -------------------------------------------------
        # evaluate p(y_t|s_t,θ), y_t = fΨ(s_t;θ) + u_t ,   u_t ∼ N(0,Σ_u(θ))

        vLogIncrWeights         = [ ll( mData[t,:], fΨ(vvs[jj]) , Σ_u )  for jj = 1:Mbspf ]


        # Likelihood increment:
        vHelp                   = vLogIncrWeights .+ log.(mWeights[t,:]) #log of w-tilde * W
        c                       = maximum(vHelp)
        vHelp                   = vHelp .- c

        vLogLik[t]              = log(mean(exp.(vHelp))) + c


        # --- c) Normalized Weights : ----------------------------------------------

        vNormWeights            = exp.(vHelp) ./ mean(exp.(vHelp))


        # --- d) Selection : --------------------------------------------------------

        ESS                     = Mbspf / mean( vNormWeights.^2 )
        doSelect                = ESS < Mbspf/2

        if doSelect == 1

            #vSelectionIndices       = rand(Categorical(vNormWeights./Mbspf),Mbspf,1)
            vSelectionIndices       = resampleNYFED(vNormWeights)
            [aParticlesEndoStates[t+1,jj,:]   = vvs[vSelectionIndices[jj]] for jj = 1:Mbspf]
            [aParticlesExoStates[t+1,jj,:]   = vvε[vSelectionIndices[jj]] for jj = 1:Mbspf]
            mWeights[t+1,:]         .= 1

        elseif doSelect == 0

            mWeights[t+1,:]         = vNormWeights
            [aParticlesEndoStates[t+1,jj,:]   = vvs[jj] for jj = 1:Mbspf]
            [aParticlesExoStates[t+1,jj,:]   = vvε[jj] for jj = 1:Mbspf]

        end

    end


    ### ------------------------ COMPUTE LL APPROX. ----------------------------

    loglik          = sum( vLogLik )


    return loglik, aParticlesEndoStates, aParticlesExoStates, mWeights, vLogLik

end



# -- Resulting Likelihood Function ---------------------------------------------


function fLL_L3(θ)

    fΦ, fΨ, fDrawInitialStatesBSPF = fSolveModelL3(θ)

    Mbspf   = 5000

    nS      = 3
    nε      = 2

    return fBSPF(mData,Mbspf,nS,nε,θ,fΨ,fΦ,fDrawInitialStatesBSPF)[1]

end


# Only used for assessing BSPF accuracy (in "main_RBC_CheckBSPF.jl"):


function fLL_L3(θ,Mbspf)

    timeSolution = @elapsed fΦ, fΨ, fDrawInitialStatesBSPF = fSolveModelL3(θ)

    nS      = 3
    nε      = 2

    timePF = @elapsed loglik = fBSPF(mData,Mbspf,nS,nε,θ,fΨ,fΦ,fDrawInitialStatesBSPF)[1]

    return loglik, timeSolution, timePF

end



# -- Function to simulate data based on this nonlinear SS rep: -----------------


function fSimulate_mL3(fΨ,fΦ,θ,nY,nS,nε, T, vs0=0, vε0=0)

    # Simulates Linear State Space of the form:
    #           s_t  =  fΦ(s_{t-1},ε_t;θ)         ,   ε_t ∼ N(0,I)
    #           y_t  =  fΨ(s_t;θ)          +  u_t ,   u_t ∼ N(0,Σ_u(θ))


    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)

    Σ_u_sqrt       = diagm([σy,σi,σl])
    nη             = 3


    mε      = zeros(nε,T+1)
    mS      = zeros(nS,T+1)
    mY      = zeros(nY,T)

    vs0 == 0 ? vs0 = [0.00,0.00,kSS] : nothing #this is [z,b,k] at SS
    vε0 == 0 ? vε0 = zeros(nε,1) : nothing

    mS[:,1] = vs0

    mu      = Σ_u_sqrt * randn(nη,T)

    for tt = 1:T
        #println(tt)
        mε[:,tt+1] = randn(2)
        mS[:,tt+1] = fΦ(mS[:,tt],mε[:,tt+1])
        mY[:,tt]   = fΨ(mS[:,tt+1]) + mu[:,tt]
    end


    return mY', mS', mε', mu'

end
