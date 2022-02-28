# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "HMpsi" (homoskedastic VAR(1) with intercept, using Φε instead of Σ)



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

vParLabs_HMpsi         = [vec(mLabs_Φ)..., mLabs_Φε[tril(trues(size(mLabs_Φε)))]...]


vGroupLengths_HMpsi    = [n*k, n*(n+1)/2]
mParaPointers_HMpsi    = fGetParaPointers(vGroupLengths_HMpsi)


# Functions to go from ϑ to VAR objects and back:

function fϑtoVARObjects_HMpsi(ϑ)

    n     = 2
    p     = 1
    k     = n*p + 1

    Φ     = reshape(ϑ[mParaPointers_HMpsi[1,1]:mParaPointers_HMpsi[1,2]],k,n)

    Φε    = zeros(n,n)
    Φε[tril(trues(size(Φε)))] = ϑ[mParaPointers_HMpsi[2,1]:mParaPointers_HMpsi[2,2]]

    return Φ, Φε

end


function fVARObjectsToϑ_HMpsi(Φ, Φε)

    ϑ   = [vec(Φ)..., Φε[tril(trues(size(Φε)))]...]

    return ϑ

end


# Log-likelihood:

ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)


function fLL_HMpsi(ϑ)

    Φ, Φε   = fϑtoVARObjects_HMpsi(ϑ)
    Σ       = Φε * Φε'
    YY, XX  = fDataVARtoReg(mData,1)
    T, n    = size(YY)

    vLL     = [ll(YY[t,:],Φ'*XX[t,:],Σ) for t=1:T]

    return sum(vLL)

end
