# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Function to simulate linear state space



# -------------------------------------------------------------------------------


function fSimulateLinearSS(Ψ0, Ψ1, Ψu, Φ1, Φε, T, X0=0)

    # Simulates Linear State Space of the form:
    #           s_t =      Φ1 s_{t-1} + Φε ε_t ,   ε_t ∼ N(0,I)
    #           y_t = Ψ0 + Ψ1 s_t     + Ψu u_t ,   u_t ∼ N(0,I)

    nY, nX  = size(Ψ1)
    nE      = size(Φε,2)

    mX      = zeros(nX,T+1)
    mY      = zeros(nY,T)

    X0 == 0 ? X0 = zeros(nX,1) : nothing

    mX[:,1] = X0

    mu      = Ψu * randn(nY,T) # drawing all the u's first is compatible with nonlinear SS, where only u's enter additively and can hence be drawn first, before the "loop" over t below. ε has to be drawn for each t separately in function...
    mε      = Φε * randn(nE,T)

    [mX[:,tt+1] =      Φ1 * mX[:,tt]   + mε[:,tt]  for tt = 1:T ]
    [mY[:,tt]   = Ψ0 + Ψ1 * mX[:,tt+1] + mu[:,tt]  for tt = 1:T ]

    return mY', mX', mε', mu'

end
