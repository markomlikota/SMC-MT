# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Prior specification file.

# relates to model "HM"

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

S_pr_HM     = deepcopy(S_pr)
v_pr_HM     = deepcopy(v_pr)

μ_pr_HM     = deepcopy(μ_pr)
P_pr_HM     = deepcopy(P_pr)



# Define fPriorDraw_Mi(), fPriorLogEval_Mi(ϑMi), fIsDrawValid_Mi(ϑMi):

function fPriorDraw_HM()

    # Draw model parameter objects:

    Σ, Φ        = fOneDrawMNIW(μ_pr_HM,P_pr_HM,S_pr_HM,v_pr_HM)


    # Assemble parameter vector ϑ:

    ϑ           = fVARObjectsToϑ_HM(Φ, Σ)


    return ϑ

end



function fPriorLogEval_HM(ϑ)

    # Extract objects from ϑ:

    Φ, Σ            = fϑtoVARObjects_HM(ϑ)


    # Compute prior probabilities of individual parameters (-groups) in ϑ:

    prprΦΣ          = fEvalMNIW(Φ,Σ,μ_pr_HM,P_pr_HM,S_pr_HM,v_pr_HM)


    # Put them together into a vector of same dimension as ϑ:

    vP              = ones(length(ϑ))

    vP[mParaPointers_HM[1,1]]           = prprΦΣ


    return log.(vP)

end



function fIsDrawValid_HM(ϑ)

    # Extract needed objects from ϑ:

    Φ, Σ           = fϑtoVARObjects_HM(ϑ)


    # Check validity for each:

    # isVARstationary =
    # try
    #     maximum(abs.(eigvals(fCompanionForm(Φ,Σ)[2]))) < 1
    # catch
    #     false
    # end

    vEigValsΣ = eigvals(Σ)
    isSigPD = ( sum(vEigValsΣ .< 0) == 0 ) && ( sum(imag.(vEigValsΣ) .!= 0) == 0 )


    # Put them together into a vector of same dimension as ϑ:

    vValid          = trues(length(ϑ))

    # vValid[mParaPointers_HM[1,1]]       = isVARstationary
    vValid[mParaPointers_HM[2,1]]         = isSigPD


    return vValid

end
