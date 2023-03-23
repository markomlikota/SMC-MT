# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Prior specification file.

# relates to  model "L1" (first order linearized RBC)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Load external functions:

# include(sMyPath * "Functions/fVAR.jl")



# Load DGPspec (because some parameters in prior functions below are fixed at their true values):

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )



# Define fPriorDraw_Mi(), fPriorLogEval_Mi(θ), fIsDrawValid_Mi(θ):


function fPriorDraw_L11()


    # Draw model parameter objects:

    ρz     = rand(Beta(fMyBeta(0.6,0.15)...))   #rand()
    σz     = rand(InverseGamma(fMyInvGamma(5.0,1.5)...))/100

    ρb     = rand(Beta(fMyBeta(0.6,0.15)...))   #rand()
    σb     = rand(InverseGamma(fMyInvGamma(5.0,1.5)...))/100

    r      = r_0
    τ      = rand(Gamma(fMyGamma(1.0,1.0)...))
    α      = rand(Beta(fMyBeta(0.35,0.05)...))
    δ      = δ_0
    ν      = ν_0

    ϕ1     = rand(Gamma(fMyGamma(30,15)...))
    ϕ2     = ϕ2_0 #rand(Normal(0,50)) #note: 2nd argument is StD, not Var!

    σy     = σy_0
    σi     = σi_0
    σl     = σl_0


    # Assemble parameter vector ϑ:

    θ           = [ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl]


    return θ

end



function fPriorLogEval_L11(θ)


    # Extract parameters from θ:

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ


    # Compute prior probabilities of individual parameters in θ:

    vP     = zeros(length(θ))


    vP[1]  = pdf(Beta(fMyBeta(0.6,0.15)...),ρz)
    vP[2]  = pdf(InverseGamma(fMyInvGamma(5.0,1.5)...),σz*100)

    vP[3]  = pdf(Beta(fMyBeta(0.6,0.15)...),ρb)
    vP[4]  = pdf(InverseGamma(fMyInvGamma(5.0,1.5)...),σb*100)

    vP[5]  = 1
    vP[6]  = pdf(Gamma(fMyGamma(1.0,1.0)...),τ)
    vP[7]  = pdf(Beta(fMyBeta(0.35,0.05)...),α)
    vP[8]  = 1
    vP[9]  = 1

    vP[10] = pdf(Gamma(fMyGamma(30,15)...),ϕ1)
    vP[11] = 1 #pdf(Normal(0,50),ϕ2)

    vP[12] = 1
    vP[13] = 1
    vP[14] = 1


    return log.(vP)

end


# fDomainCheck(x,lower,upper) = x > lower && x < upper

# vBoundTypes = [2,1,2,1,1,1,2,2,1,1,0,1,1] #needed in output analysis (Kernel Densities for PrPoPlots)

function fIsDrawValid_L11(θ)


    # Extract parameters from θ:

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ


    # Check validity for each:

    vValid = zeros(length(θ))

    vValid[1] = ρz >= 0 && ρz < 1
    vValid[2] = σz > 0

    vValid[3] = ρb >= 0 && ρb < 1
    vValid[4] = σb > 0

    vValid[5] = r > 0
    vValid[6] = τ > 0
    vValid[7] = α > 0 && α < 1
    vValid[8] = δ > 0 && δ < 1
    vValid[9] = ν > 0

    vValid[10] = ϕ1 >= 0
    vValid[11] = true

    vValid[12] = σy >= 0
    vValid[13] = σi >= 0
    vValid[14] = σl >= 0


    vValid = vValid.==1


    return vValid

end
