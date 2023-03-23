# Marko Mlikota, University of Pennsylvania, mlikota@sas.upenn.edu
# July 2020
# -------------------------------------------------------------------------------


# Functions to perform Kalman filter and compute log-likelihood.


# -------------------------------------------------------------------------------


myChol(x) = cholesky(x).L

myInv(x) = inv(x) #pinv(x; atol=1e-15)

drawRandMVN(vμ,mΣ) = rand(MvNormal(vμ,Matrix(Hermitian(mΣ))))

function fMyMNDraw(vμ,mΣ)

      n = length(vμ)
      F = svd(mΣ)

      return vμ + F.U * Diagonal(sqrt.(F.S)) * rand(MvNormal(zeros(n),I))

end


function fKalmanFilter(mData,Ψ0,Ψ1,Φ1,Φε,s0=0,P00=1,Σ_u=0)

      # Inputs: - matrix (T x nY) of data (for Y)
      #         - matrices for state space representation of the form:
      #           s_t =      Φ1 s_{t-1} + Φε ε_t ,   ε_t ∼ N(0,I)
      #           y_t = Ψ0 + Ψ1 s_t     +    u_t ,   u_t ∼ N(0,Σ_u)
      #         - initial state: s0
      #                 if = 0: taken to be zero (mean),
      #                 otherwise supply value (nS x 1 vector)
      #         - initial covariance matrix for state: P00
      #                 if = 0: matrix of zeros (no variance, initial state known)
      #                 if = 1: unconditional variance
      #         - measurement error covariance matrix Σ_u
      #                 if = 0: there are no measurement errors

      # Outputs: - 2 x 1 vector of arrays: nS x T matrix of predicted states (means), nS x nS x T array of prediction variances
      #          - analogous for updated states (2 x 1 vector with nS x (T+1) matrix and nS x nS x (T+1) array)
      #          - analogous for prediction of measurement variable y (2 x 1 vector with nS x T matrix and nS x nS x T array)
      #          - vector of log likelihood increments


      ### PRELIMINARY:
      # -------------------------------------------------

      T,nY        = size(mData)
      nS          = size(Φ1)[1]
      Σ_u == 0 ? Σ_u = zeros(nY,nY) : nothing


      ### CREATE OBJECTS TO STORE RESULTS:
      # -------------------------------------------------

      mPredMeanS  = zeros(nS,T)           # s_{t|t-1} (mean of RV s_t|Y_{1:t-1})
      aPredVarS   = zeros(nS,nS,T)        # P_{t|t-1} (variance)

      mUpdMeanS   = zeros(nS,T+1)         # s_{t|t} (mean of RV s_t|Y_{1:t})
      aUpdVarS    = zeros(nS,nS,T+1)      # P_{t|t} (variance)

      mPredMeanY  = zeros(nY,T)           # y_{t|t-1} (mean of RV y_t|Y_{1:t-1})
      aPredVarY   = zeros(nY,nY,T)        # F_{t|t-1} (variance)

      vLogLik     = zeros(T)


      ### INITIALIZE P_{0|0}, s_0:
      # -------------------------------------------------

      # Initialize P_0|0:

            # a) take to be zero matrix (known, fixed initial state):
            P00 == 0 ? P00vec = zeros(nS*nS) : nothing

            # b) unconditional variance of s_t:
            P00 == 1 ? P00vec = myInv( diagm(vec(ones(nS^2,1))) - kron(Φ1,Φ1) ) * vec(Φε * Φε') : nothing

            # c) h-step ahead conditional variance of s_t:


            # Store:

            aUpdVarS[:,:,1]   = reshape(P00vec,(nS,nS))


      # Initialize s_0 (if =0):

            s0 == 0 ? s0 = zeros(nS,1) : nothing


            # Store:

            mUpdMeanS[:,1] = s0


      ### RECURSIONS
      # -------------------------------------------------

      for t = 1:T

            # println(t)

            tt = t+1 #replaces the index "t" for objects with dimensions T+1 instead of T (m/a Upd)

            # --- 1. Forecasting ---

            ## a) States
            mPredMeanS[:,t]   = Φ1* mUpdMeanS[:,tt-1]
            aPredVarS[:,:,t]  = Φ1* aUpdVarS[:,:,tt-1]* Φ1' + Φε* Φε'

            ## b) Measurement
            mPredMeanY[:,t]   = Ψ0 + Ψ1* mPredMeanS[:,t]
            aPredVarY[:,:,t]  = Ψ1* aPredVarS[:,:,t]* Ψ1' + Σ_u

            #compute LL increment:

            vLogLik[t]        = ll(mData[t,:],mPredMeanY[:,t],aPredVarY[:,:,t])


            # --- 2. Updating ---

            predError         = mData[t,:] - mPredMeanY[:,t]
            mUpdMeanS[:,tt]   = mPredMeanS[:,t] +
                                aPredVarS[:,:,t]* Ψ1'* myInv(aPredVarY[:,:,t])* predError
            aUpdVarS[:,:,tt]  = aPredVarS[:,:,t] -
                                aPredVarS[:,:,t]* Ψ1'* myInv(aPredVarY[:,:,t])* Ψ1* aPredVarS[:,:,t]


      end

      vaPredS     = Array{Float64}[mPredMeanS, aPredVarS]
      vaUpdS      = Array{Float64}[mUpdMeanS, aUpdVarS]
      vaPredY     = Array{Float64}[mPredMeanY, aPredVarY]

      return vaPredS, vaUpdS, vaPredY, vLogLik

end

ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)

