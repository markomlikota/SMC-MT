# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Prior specification file.

# relates to model "HMpsi"

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

S_pr_HMpsi     = deepcopy(S_pr)
v_pr_HMpsi     = deepcopy(v_pr)

μ_pr_HMpsi     = deepcopy(μ_pr)
P_pr_HMpsi     = deepcopy(P_pr)



# Define fPriorDraw_Mi(), fPriorLogEval_Mi(ϑMi), fIsDrawValid_Mi(ϑMi):


mK = [1 0 0 0;
      0 0 1 0;
      0 1 0 0;
      0 0 0 1]

mD = [1 0 0;
      0 1 0;
      0 0 0;
      0 0 1]

mL = [1 0 0 0;
      0 1 0 0;
      0 0 0 1]

fJacobian(Φε) = mL * (I + mK) * kron(Φε,Matrix(I,n,n)) * mD

function fPriorDraw_HMpsi()

    # Draw model parameter objects:

    while true

        Σ, Φ        = fOneDrawMNIW(μ_pr_HMpsi,P_pr_HMpsi,S_pr_HMpsi,v_pr_HMpsi)

        Φε          = cholesky(Σ).L

        pAccept     = min( det( fJacobian(Φε) ) / 30 , 1)
        doAccept    = rand( Bernoulli(pAccept) )
        if doAccept == 1
            ϑ           = fVARObjectsToϑ_HMpsi(Φ, Φε)
            return ϑ
        end

    end

    # # Assemble parameter vector ϑ:
    #
    # ϑ           = fVARObjectsToϑ_HMpsi(Φ, Φε)


    return ϑ

end


function fPriorLogEval_HMpsi(ϑ)

    # Extract objects from ϑ:

    Φ, Φε           = fϑtoVARObjects_HMpsi(ϑ)
    Σ               = Φε * Φε'


    # Compute prior probabilities of individual parameters (-groups) in ϑ:

    prprΦΣ          = fEvalMNIW(Φ,Σ,μ_pr_HMpsi,P_pr_HMpsi,S_pr_HMpsi,v_pr_HMpsi)


    # Put them together into a vector of same dimension as ϑ:

    vP              = ones(length(ϑ))

    vP[mParaPointers_HMpsi[1,1]]           = prprΦΣ * det( fJacobian(Φε) )


    return log.(vP)

end



function fIsDrawValid_HMpsi(ϑ)

    # Extract needed objects from ϑ:

    Φ, Φε           = fϑtoVARObjects_HMpsi(ϑ)
    Σ               = Φε * Φε'


    # Check validity for each:

    # isVARstationary = maximum(abs.(eigvals(fCompanionForm(Φ,Σ)[2]))) < 1


    # Put them together into a vector of same dimension as ϑ:

    vValid          = trues(length(ϑ))

    # vValid[mParaPointers_HMpsi[1,1]]       = isVARstationary


    return vValid

end
