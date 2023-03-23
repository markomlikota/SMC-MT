# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Prior specification file.

# relates to model SV

# stationarity not imposed!



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Load external functions:

include(sMyPath * "Functions/fVAR.jl")
include(sMyPath * "Functions/fMNIW.jl")



# Construct prior MNIW parameters using Minnesota prior (dummy observations):

vλ          = [1,1,3]

y_bar       = mean(mData,dims=1)
s_bar       = std(mData,dims=1)

Ystar       = [ diagm(vec(vλ[1]*s_bar))                 ;
                vλ[2]*y_bar                             ;
                repeat(diagm(vec(s_bar)),Int(vλ[3]),1)  ]

Xstar       = [ diagm(vec(vλ[1]*s_bar)) zeros(2,1)      ;
                vλ[2]*y_bar vλ[2]                       ;
                zeros(n*Int(vλ[3]),k)                   ]

Φstar       = (Xstar'*Xstar)^(-1) * Xstar' * Ystar
Sstar       = (Ystar-Xstar*Φstar)'*(Ystar-Xstar*Φstar)
Tstar       = size(Ystar,1)


S_pr        = Sstar
v_pr        = Tstar - k

μ_pr        = Φstar
P_pr        = Xstar'*Xstar


# Hardcode prior distribution parameters for this model:

S_pr_SV     = deepcopy(S_pr)
v_pr_SV     = deepcopy(v_pr)

μ_pr_SV     = deepcopy(μ_pr)
P_pr_SV     = deepcopy(P_pr)



# Define fPriorDraw_Mi(), fPriorLogEval_Mi(ϑMi), fIsDrawValid_Mi(ϑMi):

function fPriorDraw_SV()

    n           = 2


    # Draw model parameter objects:

    Σ, Φ        = fOneDrawMNIW(μ_pr_SV,P_pr_SV,S_pr_SV,v_pr_SV)

    vρ          = zeros(n)
    [vρ[nn]     = rand(Uniform(0,1)) for nn = 1:n]

    vξ          = zeros(n)
    [vξ[nn]     = sqrt.(rand(InverseGamma(fMyInvGamma2(0.3,2.0)...))) for nn = 1:n]


    # Assemble parameter vector ϑ:

    ϑ           = fVARObjectsToϑ_SV(Φ, Σ, vρ, vξ)


    return ϑ

end


function fPriorLogEval_SV(ϑ)

    # Extract objects from ϑ:

    Φ, Σ, vρ, vξ   = fϑtoVARObjects_SV(ϑ)


    # Compute prior probabilities of individual parameters (-groups) in ϑ:

    prprΦΣ        = fEvalMNIW(Φ,Σ,μ_pr_SV,P_pr_SV,S_pr_SV,v_pr_SV)

    vprprρ        = [pdf(Uniform(0,1), vρ[nn]) for nn = 1:n]
    vprprξ        = [pdf(InverseGamma(fMyInvGamma2(0.3,2.0)...), vξ[nn]^2) for nn = 1:n]


    # Put them together into a vector of same dimension as ϑ:

    vP            = ones(length(ϑ))

    vP[mParaPointers_SV[1,1]]                       = prprΦΣ
    vP[mParaPointers_SV[3,1]:mParaPointers_SV[3,2]] = vprprρ
    vP[mParaPointers_SV[4,1]:mParaPointers_SV[4,2]] = vprprξ


    return log.(vP)

end



function fIsDrawValid_SV(ϑ)

    # Extract needed objects from ϑ:

    Φ, Σ, vρ, vξ   = fϑtoVARObjects_SV(ϑ)


    # Check validity for each:

    # isVARstationary =
    # try
    #     maximum(abs.(eigvals(fCompanionForm(Φ,Σ)[2]))) < 1
    # catch
    #     false
    # end

    vEigValsΣ = eigvals(Σ)
    isSigPD = ( sum(vEigValsΣ .< 0) == 0 ) && ( sum(imag.(vEigValsΣ) .!= 0) == 0 )

    vIsSVstationary = abs.(vρ) .< 1
    vAreXisPositive = vξ .> 0


    # Put them together into a vector of same dimension as ϑ:

    vValid          = trues(length(ϑ))
    # vValid[mParaPointers_SV[1,1]]                     = isVARstationary
    vValid[mParaPointers_SV[1,1]]                       = isSigPD
    vValid[mParaPointers_SV[3,1]:mParaPointers_SV[3,2]] = vIsSVstationary
    vValid[mParaPointers_SV[4,1]:mParaPointers_SV[4,2]] = vAreXisPositive


    return vValid

end
