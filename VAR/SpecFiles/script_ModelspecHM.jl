# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "HM" (homoskedastic VAR(1) with intercept, using Σ instead of Φε)



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

mLabs_Σ             = string.(repeat(["sig_"],n,n) .* repeat(string.(collect(1:n)),1,n) .* repeat(string.(collect(1:n)'),n,1) )

vParLabs_HM         = [vec(mLabs_Φ)..., mLabs_Σ[tril(trues(size(mLabs_Σ)))]...]


vGroupLengths_HM    = [n*k, n*(n+1)/2]
mParaPointers_HM    = fGetParaPointers(vGroupLengths_HM)


# Functions to go from ϑ to VAR objects and back:

function fϑtoVARObjects_HM(ϑ)

    n     = 2
    p     = 1
    k     = n*p + 1

    Φ     = reshape(ϑ[mParaPointers_HM[1,1]:mParaPointers_HM[1,2]],k,n)

    Σ                       = zeros(n,n)
    Σ[tril(trues(size(Σ)))] = ϑ[mParaPointers_HM[2,1]:mParaPointers_HM[2,2]]
    mIndHelp                = tril(trues(size(Σ)),-1)
    Σ[mIndHelp']            = Σ[mIndHelp]

    return Φ, Σ

end


function fVARObjectsToϑ_HM(Φ, Σ)

    ϑ   = [vec(Φ)..., Σ[tril(trues(size(Σ)))]...]

    return ϑ

end


# Log-likelihood:

ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)


function fLL_HM(ϑ)

    Φ, Σ    = fϑtoVARObjects_HM(ϑ)
    YY, XX  = fDataVARtoReg(mData,1)
    T, n    = size(YY)

    vLL     = [ll(YY[t,:],Φ'*XX[t,:],Σ) for t=1:T]

    return sum(vLL)

end
