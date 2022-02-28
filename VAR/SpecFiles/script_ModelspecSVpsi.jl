# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "SVpsi" (VAR(1) with intercept and stochastic volatility, using Φε instead of Σ)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Load external functions needed for this script:

include(sMyPath * "Functions/fVAR.jl")



# Labels and parameter-pointers under this model :

n                   = 2
p                   = 1
k                   = n*p + 1


mLabs_Φ1            = string.(repeat(["phi1_"],n,n) .* repeat(string.(collect(1:n)),1,n) .* repeat(string.(collect(1:n)'),n,1))
mLabs_Φc            = string.(repeat(["phic_"],n,1) .* string.(collect(1:n)))
mLabs_Φ             = permutedims([mLabs_Φ1 mLabs_Φc])

mLabs_Φε            = string.(repeat(["phi_eps_"],n,n) .* repeat(string.(collect(1:n)),1,n) .* repeat(string.(collect(1:n)'),n,1) )

vLabs_vρ            = string.(repeat(["rho_"],n,1) .* string.(collect(1:n)))

vLabs_vξ            = string.(repeat(["xi_"],n,1) .* string.(collect(1:n)))

vParLabs_SVpsi         = [vec(mLabs_Φ)..., mLabs_Φε[tril(trues(size(mLabs_Φε)))]..., vLabs_vρ..., vLabs_vξ...]


vGroupLengths_SVpsi    = [n*k, n*(n+1)/2, n, n]
mParaPointers_SVpsi    = fGetParaPointers(vGroupLengths_SVpsi)


# Functions to go from ϑ to VAR objects and back:

function fϑtoVARObjects_SVpsi(ϑ)

    n     = 2
    p     = 1
    k     = n*p + 1

    Φ     = reshape(ϑ[mParaPointers_SVpsi[1,1]:mParaPointers_SVpsi[1,2]],k,n)

    Φε    = zeros(n,n)
    Φε[tril(trues(size(Φε)))] = ϑ[mParaPointers_SVpsi[2,1]:mParaPointers_SVpsi[2,2]]

    vρ    = ϑ[mParaPointers_SVpsi[3,1]:mParaPointers_SVpsi[3,2]]

    vξ    = ϑ[mParaPointers_SVpsi[4,1]:mParaPointers_SVpsi[4,2]]

    return Φ, Φε, vρ, vξ

end


function fVARObjectsToϑ_SVpsi(Φ, Φε, vρ, vξ)

    ϑ   = [vec(Φ)..., Φε[tril(trues(size(Φε)))]..., vρ..., vξ...]

    return ϑ

end



# Bootstrap particle filter to compute LL:

function fBSPF(ϑ,mData,Mbspf)

    # Bootstrap Particle Filter

    # Inputs:
    # - vector of parameters (to be used in proposal, likelihood and initialization densities)
    # - (T x nY) matrix of data for likelihood evaluation,
    # - number of particles

    # Output: loglik, array of particles, matrix of weights



    ### ------------------------ SET UP ----------------------------------------

    Φ, Φε, vρ, vξ   = fϑtoVARObjects_SVpsi(ϑ)

    ns              = length(vρ)
    DD              = Matrix(I,ns,ns)

    #minimum(vξ .== 0) == 1 ? Mbspf = 1 : nothing # for homoskedastic model (VAR); to test whether function works


    YY, XX          = fDataVARtoReg(mData,1)
    T, n            = size(YY)



    ### ------------------------ INITIALIZE OBJECTS ----------------------------

    aParticles      = zeros(T+1,Mbspf,ns)
    mWeights        = zeros(T+1,Mbspf)
    vLogLikIncr     = zeros(T,1)



    ### ------------------------ INITIALIZE ALGORITHM (RECURSIONS) -------------


    # Take one-step ahead density to draw s0:
    # p(s_t|s_{t-1}=DD,ϑ)

    #vvs0                 = [fDrawProposal(diag(log.(DD)))' for m = 1:M ]


    # Take unconditional density, or h-step ahead density in case of unit root:
    # p(s_t|ϑ), or p(s_t|s_{t-h}=DD,ϑ), h=20

    ρmax                 = maximum(vρ)
    if ρmax < 0.98
        vvs0             = [fDrawUncond(DD,vρ,vξ)' for m = 1:Mbspf]
    else
        vvs0             = [fDrawSimForward(diag(log.(DD)),20,DD,vρ,vξ)' for m = 1:Mbspf ]
    end


    [aParticles[1,m,:]   = vvs0[m] for m = 1:Mbspf]
    mWeights[1,:]       .= 1



    ### ------------------------ RECURSIONS ------------------------------------

    for t = 1:T


        # --- a) Forecasting s_t : -------------------------------------------------

        vvs_                     = [ fDrawProposal(aParticles[t,jj,:],DD,vρ,vξ) for jj = 1:Mbspf ]


        # --- b) Forecasting y_t : -------------------------------------------------

        vYY                     = YY[t,:]
        vXX                     = XX[t,:]
        vLogIncrWeights         = [fEvalLogLikIncr(vvs_[jj],vYY,vXX,Φ,Φε) for jj = 1:Mbspf ]

        # Likelihood increment:
        vHelp                   = vLogIncrWeights .+ log.(mWeights[t,:]) #log of w-tilde * W
        c                       = maximum(vHelp)
        vHelp                   = vHelp .- c

        vLogLikIncr[t]          = log(mean(exp.(vHelp))) + c


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

    loglik          = sum( vLogLikIncr )


    return loglik, aParticles, mWeights, vLogLikIncr

end


function fDrawProposal(vs,DD,vρ,vξ)

    # p(s_t|s_{t-1},ϑ), input is s_{t-1}

    vMean       = log.(diag(DD)) + vρ .* (vs - log.(diag(DD)) )
    mVar        = diagm(vξ.^2)

    if minimum(vξ .== 0) == 1
        return vMean
    else
        return rand(MvNormal(vMean,mVar))
    end

end

function fDrawUncond(DD,vρ,vξ)

    # p(s_t|ϑ)

    vMean       = log.(diag(DD))
    mVar        = diagm(vξ.^2 ./(1 .- vρ))

    if minimum(vξ .== 0) == 1
        return vMean
    else
        return rand(MvNormal(vMean,mVar))
    end

end

function fDrawSimForward(vs,h,DD,vρ,vξ)

    # p(s_t|s_{t-h},ϑ), inputs are s_{t-h} & h

    vMean       = log.(diag(DD)) + vρ .* (vs .- log.(diag(DD)))
    mVar        = diagm(vξ.^2 .* (1 .- vρ.^(2*h))./(1 .- vρ))

    if minimum(vξ .== 0) == 1
        return vMean
    else
        return rand(MvNormal(vMean,mVar))
    end

end


ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)

# vs_=  vvs_[jj]
function fEvalLogLikIncr(vs_,vYY,vXX,Φ,Φε)

    # log( p(y_t|s_t,Y_{1:t-1},ϑ) ), input is s_t, and data for y_t & x_t

    DDt              = diagm(exp.(vs_))
    Σt               = Φε * DDt * Φε'

    return ll(vYY, Φ' * vXX, Σt)

end


# Final function that computes LL under this model:

fLL_SVpsi(ϑ) = fBSPF(ϑ,mData,100)[1]
