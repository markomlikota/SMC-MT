# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# July 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Prior specification file.

# relates to  model "L2" (it's the same for any non-linear model solution (i.e. cases where adjustment cost asymmetry parameters are identified))



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Load external functions:

include(sMyPath * "/Functions/fHelpersSMC.jl")



# Load DGPspec (because some parameters in prior functions below are fixed at their true values):

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )



# Define fPriorDraw_Mi(), fPriorLogEval_Mi(θ), fIsDrawValid_Mi(θ):


function fPriorDraw_L2()

    rA      = rand(Gamma(fMyGamma(2.0,1.0)...)) # defines β     = 1/(1+rA/400)
    πA      = rand(Gamma(fMyGamma(3.0,1.0)...)) # defines πStar = 1 + πA/400
    γA      = rand(Gamma(fMyGamma(2.0,1.5)...)) # defines γ     = 1 + γA/400

    τ       = rand(Gamma(fMyGamma(2.0,1.0)...))
    ν       = rand(Gamma(fMyGamma(0.5,1.0)...))

    κ       = rand(Gamma(fMyGamma(0.3,0.2)...)) # defines φp    = (1-λpSS)/(λpSS*κ*πStar^2)
    φw      = rand(Gamma(fMyGamma(15.0,7.5)...))
    ψw      = rand(Uniform(-200,200))
    ψp      = rand(Uniform(-300,300))

    ψ1      = rand(Gamma(fMyGamma(1.5,0.5)...))
    ψ2      = rand(Gamma(fMyGamma(0.2,0.1)...))

    ρr      = rand(Beta(fMyBeta(0.5,0.2)...))
    ρg      = rand(Beta(fMyBeta(0.8,0.1)...))
    ρa      = rand(Beta(fMyBeta(0.2,0.1)...))
    ρp      = rand(Beta(fMyBeta(0.6,0.2)...))

    σr      = rand(InverseGamma(fMyInvGamma(0.2,2.0)...))/100
    σg      = rand(InverseGamma(fMyInvGamma(0.75,2.0)...))/100
    σa      = rand(InverseGamma(fMyInvGamma(0.75,2.0)...))/100
    σp      = rand(InverseGamma(fMyInvGamma(0.75,2.0)...))/100

    gStar   = gStar_0
    λpSS    = λpSS_0
    λwSS    = λwSS_0
    χh      = χh_0

    σ_dy    = σ_dy_0
    σ_dwnom = σ_dwnom_0
    σ_r     = σ_r_0
    σ_infl  = σ_infl_0


    # Assemble parameter vector ϑ:

    θ          = [rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, σ_dy, σ_dwnom, σ_r, σ_infl]


    return θ

end


function fPriorLogEval_L2(θ)


    # Extract parameters from θ:

    rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, σ_dy, σ_dwnom, σ_r, σ_infl = θ


    # Compute prior probabilities of individual parameters in θ:

    vP     = ones(length(θ))

    vP[1]      = pdf(Gamma(fMyGamma(2.0,1.0)...),rA) # defines β     = 1/(1+rA/400)
    vP[2]      = pdf(Gamma(fMyGamma(3.0,1.0)...),πA) # defines πStar = 1 + πA/400
    vP[3]      = pdf(Gamma(fMyGamma(2.0,1.5)...),γA) # defines γ     = 1 + γA/400

    vP[4]       = pdf(Gamma(fMyGamma(2.0,1.0)...),τ)
    vP[5]       = pdf(Gamma(fMyGamma(0.5,1.0)...),ν)

    vP[6]       = pdf(Gamma(fMyGamma(0.3,0.2)...),κ) # defines φp    = (1-λpSS)/(λpSS*κ*πStar^2)
    vP[7]      = pdf(Gamma(fMyGamma(15.0,7.5)...),φw)
    vP[8]      = pdf(Uniform(-200,200),ψw)
    vP[9]      = pdf(Uniform(-300,300),ψp)

    vP[10]      = pdf(Gamma(fMyGamma(1.5,0.5)...),ψ1)
    vP[11]      = pdf(Gamma(fMyGamma(0.2,0.1)...),ψ2)

    vP[12]      = pdf(Beta(fMyBeta(0.5,0.2)...),ρr)
    vP[13]      = pdf(Beta(fMyBeta(0.8,0.1)...),ρg)
    vP[14]      = pdf(Beta(fMyBeta(0.2,0.1)...),ρa)
    vP[15]      = pdf(Beta(fMyBeta(0.6,0.2)...),ρp)

    vP[16]      = pdf(InverseGamma(fMyInvGamma(0.2,2.0)...),σr*100)
    vP[17]      = pdf(InverseGamma(fMyInvGamma(0.75,2.0)...),σg*100)
    vP[18]      = pdf(InverseGamma(fMyInvGamma(0.75,2.0)...),σa*100)
    vP[19]      = pdf(InverseGamma(fMyInvGamma(0.75,2.0)...),σp*100)

    vP[20]      = 1
    vP[21]      = 1
    vP[22]      = 1
    vP[23]      = 1

    vP[24]      = 1
    vP[25]      = 1
    vP[26]      = 1
    vP[27]      = 1

    return log.(vP)

end


# fDomainCheck(x,lower,upper) = x > lower && x < upper

# vBoundTypes = [2,1,2,1,1,1,2,2,1,1,0,1,1] #needed in output analysis (Kernel Densities for PrPoPlots)

function fIsDrawValid_L2(θ)


    # Extract parameters from θ:

    rA, πA, γA, τ, ν, κ, φw, ψw, ψp, ψ1, ψ2, ρr, ρg, ρa, ρp, σr, σg, σa, σp, gStar, λpSS, λwSS, χh, σ_dy, σ_dwnom, σ_r, σ_infl = θ


    # Check validity for each:

    vValid = ones(length(θ))

    vValid[1] = rA > 0
    vValid[2] = πA > 0
    vValid[3] = γA > 0

    vValid[4] = τ > 0
    vValid[5] = ν > 0

    vValid[6] = κ > 0
    vValid[7] = φw > 0
    vValid[8] = ψw > -200 && ψw < 200
    vValid[9] = ψp > -200 && ψp < 200

    vValid[10] = ψ1 > 0
    vValid[11] = ψ2 > 0

    vValid[12] = ρr > 0 && ρr < 1
    vValid[13] = ρg > 0 && ρg < 1
    vValid[14] = ρa > 0 && ρa < 1
    vValid[15] = ρp > 0 && ρp < 1

    vValid[16] = σr > 0
    vValid[17] = σg > 0
    vValid[18] = σa > 0
    vValid[19] = σp > 0


    vValid = vValid.==1


    return vValid

end
