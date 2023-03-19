# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Bootstrap Particle Filter (BSPF) for linear state space model



# -------------------------------------------------------------------------------


myInv(x) = inv(x) #pinv(x; atol=1e-15)

ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)

function fMyMNDraw(vμ,mΣ)

      n = length(vμ)
      F = svd(mΣ)

      return vμ + F.U * Diagonal(sqrt.(F.S)) * rand(MvNormal(zeros(n),I))

end


# mData,Mbspf,Ψ0,Ψ1,Φ1,Φε,s0,P00,Σ_u = mData,Mbspf,Ψ0,Ψ1,Φ1,Φε,1,1,Ψu*Ψu'

function fBSPF_linearSS(mData,Mbspf,Ψ0,Ψ1,Φ1,Φε,P00=1,Σ_u=0)

    # Bootstrap Particle Filter for linear state space of the form
    #           s_t =      Φ1 s_{t-1} + Φε ε_t ,   ε_t ∼ N(0,I)
    #           y_t = Ψ0 + Ψ1 s_t     +    u_t ,   u_t ∼ N(0,Σ_u)
    #
    # Inputs:
    # - (T x nY) matrix of data for likelihood evaluation,
    # - number of particles
    # - SS objects
    # - initial covariance matrix for state: P00
    #         if = 1: unconditional variance
    #          (other options not supported yet (e.g. h-step ahead cond. variance))
    # - measurement error covariance matrix Σ_u
    #         if = 0: there are no measurement errors
    #
    # Output: loglik, array of particles, matrix of weights
    #


    ### ------------------------ SET UP ----------------------------------------

    T,nY        = size(mData)
    nS          = size(Φ1)[1]
    Σ_u == 0 ? Σ_u = zeros(nY,nY) : nothing


    ### ------------------------ INITIALIZE OBJECTS ----------------------------

    aParticles      = zeros(T+1,Mbspf,nS)
    mWeights        = zeros(T+1,Mbspf)
    vLogLik         = zeros(T)



    ### ------------------------ INITIALIZE ALGORITHM (RECURSIONS) -------------


    # Initialize P_0|0:

          # a) unconditional variance of s_t:
          P00 == 1 ? P00vec = myInv( diagm(vec(ones(nS^2,1))) - kron(Φ1,Φ1) ) * vec(Φε * Φε') : nothing

          # b) h-step ahead conditional variance of s_t:


          # Store:

          P00               = reshape(P00vec,(nS,nS))


    # Initialize s_0 (if =0):

          vvs0             = [fMyMNDraw(zeros(nS),P00) for jj = 1:Mbspf]


    # Store:

    [aParticles[1,jj,:]   = vvs0[jj] for jj = 1:Mbspf]
    mWeights[1,:]       .= 1



    ### ------------------------ RECURSIONS ------------------------------------

    for t = 1:T


        # --- a) Forecasting s_t : -------------------------------------------------
        # draw from p(s_t|s_{t-1},θ)

        vvs_                     = [ fMyMNDraw( Φ1* aParticles[t,jj,:] , Φε * Φε' ) for jj = 1:Mbspf ]


        # --- b) Forecasting y_t : -------------------------------------------------
        # evaluate p(y_t|s_t,θ)

        vLogIncrWeights         = [ ll( mData[t,:], vec( Ψ0 + Ψ1 * vvs_[jj] ) , Σ_u )  for jj = 1:Mbspf ]


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
            [aParticles[t+1,jj,:]   = vvs_[vSelectionIndices[jj]] for jj = 1:Mbspf]
            mWeights[t+1,:]         .= 1

        elseif doSelect == 0

            mWeights[t+1,:]         = vNormWeights
            [ aParticles[t+1,jj,:]  = vvs_[jj] for jj = 1:Mbspf ]

        end

    end


    ### ------------------------ COMPUTE LL APPROX. ----------------------------

    loglik          = sum( vLogLik )


    return loglik, aParticles, mWeights, vLogLik

end
