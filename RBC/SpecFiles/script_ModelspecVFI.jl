# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Model specification file:

# Model "RBCL1" (RBC + linex cap. adj costs, first order linearization around SS),
# where LL is evaluated using BSPF



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# -- Preliminary ---------------------------------------------------------------


# Load external functions needed for this script:

include(sMyPath * "Functions/fModelFunctions.jl")
include(sMyPath * "Functions/fTauchen.jl")
include(sMyPath * "Functions/fMyMNDraw.jl")


# Parameter labels (under this model) :

vParLabs    = ["rho_z","sig_z","rho_b","sig_b","r","tau","alpha","delta","nu","phi1","phi2","sigma_y","sigma_i","sigma_l"]



# -- Functional Forms ----------------------------------------------------------


# Utility:

function fU(B,C,L,τ,ν)

    τ == 1 && return log(C) - B * L^(1+1/ν) / (1+1/ν)

    return (C^(1-τ) - 1)/(1-τ) - B * L^(1+1/ν) / (1+1/ν)

end


# Production Function:

fY(Z,K,L,α)     = Z * K^α * L^(1-α)


# Capital Adjustment Costs:

function fAC(X,ϕ1,ϕ2)

    ϕ2 != 0 && return ϕ1 / ϕ2^2 * ( exp( -ϕ2 * (X-1) ) + ϕ2 * (X-1) - 1 )
    ϕ2 == 0 && return ϕ1 / 2 * (X-1)^2

end


# Investment:

fI(K_,K,δ) = K_ - (1-δ) * K


# Optimal Labor Supply Choice :

function fLopt(B,Z,C,K,α,ν,τ)

    term = (1-α)/B * C^(-τ) * Z * K^α

    return term^(ν/(α*ν+1))

end


# Period Return:

punishUtility = -1e10

function fGoodGuessLopt(K,K_,vs,θ)

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    Z, B = vs

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)


    I = fI(K_,K,δ)

    z    = log(Z/Zstar)
    b    = log(B/Bstar)
    khat = K/kSS - 1
    ihat = I/iSS - 1

    chat = ( cSS/ySS + τ*ν*(1-α)/(α*ν+1) )^(-1) * ( (ν+1)/(α*ν+1)*(z+α*khat) - (1-α)*ν/(α*ν+1)*b - iSS/ySS*ihat )

    lhat = ν/(α*ν+1) * ( z - b + α*khat - τ*chat )

    return (1+lhat)*lSS

end

function fPeriodReturn(K,K_,vs,θ)


    # Extract parameters from θ:

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    # Extract exogenous states:

    Z, B = vs


    # Given choice of K', compute C, L :

    fC(L)   = fY(Z,K,L,α) - fI(K_,K,δ) - K * fAC(K_/K,ϕ1,ϕ2)

    function fSolveForL(L)

        L < 0 && return 1e10

        C = fC(L)

        C < 0 && return 1e10

        return L - fLopt(B,Z,C,K,α,ν,τ)

    end

    # # CHECKS:
    # pp = plot(0.01:0.01:5,fSolveForL)
    # ylims!((-10,10))
    # display(pp)
    # θ = ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl

    # @time L  = find_zero(fSolveForL,lSS) # option: ,atol=0.01
    # @time L  = optimize(L->abs(fSolveForL(L)),lSS*0.5,lSS*1.5,GoldenSection())

    # nlsolve(fSolveForL,[1.0])

    lSS = fGetSteadyState(θ)[3]

    guessL = fGoodGuessLopt(K,K_,vs,θ)
    guessL < 0 || abs(fSolveForL(guessL))>abs(fSolveForL(lSS))  ? guessL = lSS : nothing
    abs(fSolveForL(guessL)) == 1e10 ? guessL *= 2 : nothing # large negative ϕ2s move region of L-space where solution exists (C>0) to the right...
    abs(fSolveForL(guessL)) == 1e10 ? guessL *= 2 : nothing
    abs(fSolveForL(guessL)) == 1e10 ? guessL *= 2 : nothing
    abs(fSolveForL(guessL)) == 1e10 ? guessL *= 2 : nothing
    abs(fSolveForL(guessL)) == 1e10 ? guessL *= 2 : nothing

    L  = find_zero(fSolveForL,guessL) # option: ,atol=0.01
    C  = fC(L)


    # Check domains and return utility:

    C < 0 && return punishUtility
    L < 0 && return punishUtility

    return fU(B,C,L,τ,ν)

end



# -- Functions For Solution Method ---------------------------------------------


# Function to get discretized exogenous state variable(s) (s=(Z,B)):

function fGetCartProd(vSets,vIndices=0)

    # Computes Cartesian product of the sets in vSets with indices in vIndices

    # Inputs: vector of one-dimensional arrays containing the sets, optional: vector of indices (one-dimensional Array) (default: all sets)
    # Output: vector of tuples (all possible points in resulting space)

    if vIndices==0
        sets        = vSets
    else
        sets = [vSets[vIndices[i]] for i in 1:length(vIndices)]
    end

    tSpace      = unique(collect( Iterators.product(sets...) )) #tuple
    mSpace      = [tSpace[ind][i] for ind = 1:length(tSpace), i = 1:length(sets)] #matrix

    return mSpace, tSpace

end

function fGetvvs(nz,Zstar,ρz,σz,nb,Bstar,ρb,σb)


    # Discretize processes:

    vz, mPz              = fTauchen(nz,0,ρz,σz)

    vZ  = Zstar * exp.(vz)


    vb, mPb              = fTauchen(nb,0,ρb,σb)

    vB  = Bstar * exp.(vb)


    # Combine:

    vs  = fGetCartProd([vZ,vB])[1]
    mPs = kron(mPz,mPb)
    ns  = size(vs,1)

    vvs = [vs[ss,:] for ss=1:ns]

    return vvs, mPs

end


# Function to get grid for endogenous state variable (K):

fGetKgrid(kSS,nk,fracEitherSide) = collect(range(kSS*(1-fracEitherSide),kSS*(1+fracEitherSide),length=nk))
fGetKgridAsRange(kSS,nk,fracEitherSide) = range(kSS*(1-fracEitherSide),kSS*(1+fracEitherSide),length=nk)


# Functions for Value Function Iteration (VFI):

function VFI(vkgrid,v_init,vvs,mPs,θ,sMethod,maxIt,tolerance,applyMax=1,messages=0)

    # Performs value function iteration (VFI)

    # Inputs:
    # - grid for k (vector)
    # - initial guess for v (matrix; nk x ns)
    # - method is either "GridSearch" (searching for optimal k' over the same grid for k(')) or "Maximization" (actually maximizing function at each state (gridpoint))
    # - accelerator (how often max operator should be applied; e.g. 10 means every 10th iteration)

    # Output: value function, policy funcion (indices) for k'

    #Basic structure taken from:
    #https://github.com/jesusfv/Comparison-Programming-Languages-Economics/blob/master/RBC_Julia.jl


    # Preliminary:

    # Extract β from θ :
    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    β               = 1/(1 + r/400)


    # Initialize objects:

    nk, ns          = size(v_init)

    v               = deepcopy(v_init) #value function
    #(deepcopy since ow the externally supplied v_init changes to v as we run this function)
    v_              = zeros(nk,ns) #value function next period
    gk              = zeros(nk,ns) #policy function for k
    e_v_            = zeros(nk,ns) #expected value function next period

    sMethod=="Maximization" ? itp_k = interpolate(vkgrid,BSpline(Linear())) : nothing


    # Convergence Settings:

    maxDifference   = 10.0
    # maxIt         = 100000
    # tolerance     = 10e-6
    iteration       = 0


    # Main Iterations:

    # println(ω)
    # println("")
    # println(fL(1,1,2,α,ω))
    # println("")
    #
    # println(fPeriodReturnWithFixedCost(2.3,2.4,[1,1],θ))
    # println(fPeriodReturnNoInvest(2,[1,1],θ))
    # println(fPeriodReturn(2.3,2.4,[1,1],θ))
    #throw("done")

    while(maxDifference > tolerance && iteration<maxIt)

        e_v_ = v*mPs'
        #println(gk[50,17:19]) #println(e_v_[1,1:2])
        # println(mPs[1,1:2])
        sMethod=="Maximization" ? itp_e_v_ = interpolate(e_v_, BSpline(Linear())) : nothing

        for ss in 1:ns
            vs = vvs[ss]


            # We start from previous choice (monotonicity of policy function)
            gridCapitalNextPeriod = 1


            for kk in 1:nk
                k = vkgrid[kk]

                if sMethod=="GridSearch"
                    if mod(iteration,applyMax)==0 #|| iteration <= 50 #iteration==0 #if we apply max

                        valueHighSoFar = -Inf#-1e20
                        capitalChoice  = 1

                        for kk_ in gridCapitalNextPeriod:nk
                                                                 #println(string("ss=",ss,"kk=",kk,"kk_=",kk_))
                            k_ = vkgrid[kk_]

                            # periodreturn = fPeriodReturnWithFixedCost(k,k_,vs,θ)
                            periodreturn = fPeriodReturn(k,k_,vs,θ)

                            valueProvisional = periodreturn + β*e_v_[kk_,ss]

                            if (valueProvisional>valueHighSoFar)
                        	       valueHighSoFar = valueProvisional
                        	       capitalChoice = kk_
                        	       gridCapitalNextPeriod = kk_
                            else
                        	       break # We break when we have achieved the max, concavity of value function
                            end

                        end
                        # if kk==1 && iteration == 0 && ss==1
                        #     println(fPeriodReturnWithFixedCost(k,k*1.01,[1,1]))
                        #     println(fPeriodReturnNoInvest(k,[1,1]))
                        #     println(fPeriodReturn(k,k*1.01,[1,1]))
                        # end

                        # periodreturnNoInvest = fPeriodReturnNoInvest(k,vs,θ)
                        #
                        # valueNoInvest = periodreturnNoInvest + β*e_v_[kk,ss]

                        # if valueHighSoFar > valueNoInvest
                            v_[kk,ss] = valueHighSoFar
                            gk[kk,ss] = capitalChoice
                        # else
                        #     v_[kk,ss] = valueNoInvest
                        #     gk[kk,ss] = kk
                        # end

                        # v_[kk,ss] = max(valueHighSoFar,valueNoInvest)
                        # valueHighSoFar > valueNoInvest ? gk[kk,ss] = capitalChoice : gk[kk,ss] = kk

                    else

                        kk_ = Int(gk[kk,ss])
                        v_[kk,ss] = fPeriodReturn(k,vkgrid[kk_],vs,θ) + β*e_v_[kk_,ss]

                    end
                elseif sMethod=="Maximization"
                    if mod(iteration,applyMax)==0 #|| iteration <= 50 #iteration==0

                        thisperiodreturn(k_) = -1*( periodreturn = fPeriodReturn(k,k_,vs,θ) + β*itp_e_v_(k_,ss) )

                        lower = gridCapitalNextPeriod
                        upper = nk
                        res = optimize(thisperiodreturn,lower,upper,GoldenSection()) #golden section faster than brent here
                        kk_ = Optim.minimizer(res)

                        # optimize(thisperiodreturn,[kk])

                        valueInvest = -Optim.minimum(res)
                        gridCapitalNextPeriod = kk_


                        # periodreturnNoInvest = fPeriodReturnNoInvest(k,vs,θ)
                        #
                        # valueNoInvest = periodreturnNoInvest + β*itp_e_v_(kk,ss)
                        #
                        # if valueInvest > valueNoInvest
                            v_[kk,ss] = valueInvest
                            gk[kk,ss] = kk_
                        # else
                        #     v_[kk,ss] = valueNoInvest
                        #     gk[kk,ss] = kk
                        # end

                    else

                        kk_ = gk[kk,ss]
                        v_[kk,ss] = fPeriodReturn(k,itp_k(kk_),vs,θ) + β*itp_e_v_(kk_,ss)

                    end
                else
                    throw(DomainError(sMethod, "not an acceptable method"))
                end

            end

        end

        maxDifference = maximum(abs.(v_-v))
        v, v_ = v_, v #changing v_, too, since ow maxDifference above will change to zero, and this is more efficient than v = copy(v_)

        iteration = iteration+1
        if messages==1 && (mod(iteration,10)==0 || iteration==1)
            println(" Iteration = ", iteration, " Sup Diff = ", maxDifference)
        end

    end

    messages == 1 ? println(" Iteration = ", iteration, " Sup Diff = ", maxDifference) : nothing


    return v, gk, maxDifference, iteration

end

function multigridVFI(vkgridInit,v_init,vvs,mPs,θ,nkvec,sMethod,vMaxIt,vTolerance,applyMax=1,messages=0)

    # Performs VFI iteratively, length(nkvec)+1-times, with increasing grid lengths for end. state

    # Inputs:
    # - initial grid for k (vector)
    # - initial guess for v (matrix)
    # - vector specifying number of gridpoints to use in grid-updates
    # - method is either "GridSearch" or "Maximization"
    # - optional accelerator parameter
    # - optional message-parameter (for messages of underlying VFIs)

    # Output:
    # - value function (matrix)
    # - policy funcion (matrix, indices) for k'
    # - final grid for k

    v               = deepcopy(v_init)
    vkgrid          = deepcopy(vkgridInit)
    nk, ns          = size(v_init)

    v_              = zeros(nk,ns)
    nk_             = deepcopy(nk)
    vkgrid_         = deepcopy(vkgrid)

    kmin            = vkgridInit[1]
    kmax            = vkgridInit[end]

    # percent1 = 20
    # percent2 = 50
    # prop1    = 0.5
    # prop2    = 0.3
    # prop3    = 0.2

    interps  = 0


    #nkvec = [nkvec; 0]

    while true

        interps += 1

        if messages == 1
            println(" ------------------------------------------------ ")
            println(" --- VFI ", interps, ", nk = ", nk)
            println(" ------------------------------------------------ ")
        end

        #Perform iterations:
        #interps == 1 ? nothing : v_init = vNew

        vSol, gk, maxDifference, iteration = VFI(vkgrid,v,vvs,mPs,θ,sMethod,vMaxIt[interps],vTolerance[interps],applyMax,messages)

        nk == nkvec[end] ? (return vSol, gk, vkgrid, maxDifference, iteration) : nothing

        itp_v = LinearInterpolation((vkgrid,1:ns), vSol, extrapolation_bc=Line())

        #Update grid, and interpolate new initial guess:
        nk_ = nkvec[interps]
        # nk1     = Int(prop1*nk_)
        # nk2     = Int(prop2*nk_)
        # nk3     = Int(prop3*nk_)
        # vkgrid_ = fGetKGrid(kmin,kSS,nk1,nk2,nk3,percent1,percent2)
        vkgrid_ = collect(range(kmin,kmax,length=nk_))

        v_    = map(itp_v,repeat(vkgrid_,1,ns),repeat(collect(1:ns)',nk_,1))

        nk,nk_ = nk_, nk
        vkgrid,vkgrid_ = vkgrid_,vkgrid
        v, v_   = v_, v


    end

end


# Function that solves model at given θ (i.e. gets grid for states, performs VFI and casts solution (decision rule for K') as function)


function fValueFunctionL1(K,Z,B,θ)

    # Extract parameters and SS values from θ:

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    β               = 1/(1 + r/400)

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)

    z = log(Z/Zstar)
    b = log(B/Bstar)


    return fU(Bstar,cSS,lSS,τ,ν)/(1-β) + cSS^(-τ) / β * (K-kSS) + cSS^(-τ) * ySS / (1-β*ρz) * z - lSS^(1+1/ν) / (1+1/ν) / (1-β*ρb) * b

end

function fGoodGuessValueFunction(vkgridInit,vvs,θ)

    vZ = first.(vvs)
    vB = getindex.(vvs,2)

    mK = repeat(vkgridInit,outer=(1,length(vB)))
    mZ = repeat(vZ',outer=(length(vkgridInit),1))
    mB = repeat(vB',outer=(length(vkgridInit),1))

    v_init = map((K,Z,B,)->fValueFunctionL1(K,Z,B,θ),mK,mZ,mB)

    return v_init

    #v_init = repeat(vkgridInit,outer=(1,ns))#zeros(nk,ns) .+ kSS

end

function fSolveModelVFI(θ,returnInfo=0)


    # Extract parameters and SS values from θ:

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)


    # Discretize exogenous states:

    nz = 5
    nb = 5

    vvs, mPs = fGetvvs(nz,Zstar,ρz,σz,nb,Bstar,ρb,σb)

    ns = length(vvs)


    # Compute discretized solution via VFI:

    nk              = 20 #initial grid length for capital
    fracEitherSide  = 0.2
    vkgridInit      = fGetKgrid(kSS,nk,fracEitherSide)
    v_init          = fGoodGuessValueFunction(vkgridInit,vvs,θ)

    nkvec           = [50,100,200,1000]
    vMaxIt          = 1e6*[1,1,1,1,1]
    vTolerance      = [1e-2,1e-2,1e-2,1e-2,1e-2]

    vSol, gkSol, vkgridSol, maxDifference, iteration = multigridVFI(vkgridInit,v_init,vvs,mPs,θ,nkvec,"GridSearch",vMaxIt,vTolerance,30)


    # Create interpolation object based on policy function:

    gk = map(x->vkgridSol[Int(x)],gkSol)
    nk = length(vkgridSol)
    mk_ = reshape(gk,nk,nz,nb)

    itp_k_ = extrapolate( interpolate(mk_, BSpline(Cubic(Line(OnGrid()))) ), Line())
    # @time itp_k_ = interpolate(mk_, BSpline(Cubic(Line(OnGrid()))))
    # @time itp_k_ = extrapolate( interpolate(mk_, BSpline(Cubic(Line(OnGrid()))) ), Line())

    vZ = unique(first.(vvs))
    vB = unique(getindex.(vvs,2))
    itp_kgrid = LinearInterpolation(vkgridSol, 1:nk, extrapolation_bc = Line())
    itp_Zgrid = LinearInterpolation(vZ, 1:nz, extrapolation_bc = Line())
    itp_Bgrid = LinearInterpolation(vB, 1:nb, extrapolation_bc = Line())

    function fSol_K_(K,Z,B)

        Kind = itp_kgrid(K)
        Zind = itp_Zgrid(Z)
        Bind = itp_Bgrid(B)

        return itp_k_(Kind,Zind,Bind)

    end

    if returnInfo == 1
        return fSol_K_, maxDifference, iteration
    end

    return fSol_K_

end



# -- Functions To Evaluate LL given model solution -----------------------------


# Objects for nonlinear state space representation:

function fFε(_vε,ρz,σz,ρb,σb)

    # Propagates forward exogenous states: _vε to vε

    _z, _b  = _vε

    z       = rand(Normal(ρz*_z,σz))
    b       = rand(Normal(ρb*_b,σb))

    vε      = [z,b]

    return vε

end

# _vs,vε = aParticlesEndoStates[t,jj,:],vvε[jj]
function fΦ(_vs,vε,fSol_K_,θ,Zstar,Bstar)

    # Propagates forward endogenous states, given exogenous states vε: _vs to vs



    # Extract parameter values from θ:

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)


    # Extract yesterday's endogenous states _vs and today's exogenous states vε:

    K, _L, _Y, _I   = _vs

    z, b            = vε

    Z               = Zstar * exp(z)
    B               = Bstar * exp(b)


    # Compute today's endogenous states:

    K_              = fSol_K_(K,Z,B)

    # Given choice of K', compute L :

    fC(L)           = fY(Z,K,L,α) - fI(K_,K,δ) - K * fAC(K_/K,ϕ1,ϕ2)

    function fSolveForL(L)

        L < 0 && return 1e10

        C = fC(L)

        C < 0 && return 1e10

        return L - fLopt(B,Z,C,K,α,ν,τ)

    end

    guessL = fGoodGuessLopt(K,K_,[Z,B],θ)
    guessL < 0 || abs(fSolveForL(guessL))>abs(fSolveForL(lSS))  ? guessL = lSS : nothing

    L           = find_zero(fSolveForL,guessL) # option: ,atol=0.01

    Y           = fY(Z,K,L,α)
    I           = fI(K_,K,δ)


    vs          = [K_, L, Y, I]

    return vs

end

function fΨ(vs)

    # Computes measurement variables given vector of (endogenous) state variables.

    K_, L, Y, I = vs

    vy          = log.( [Y,I,L])

    return vy

end


# BSPF:

myInv(x) = inv(x) #pinv(x; atol=1e-15)

ll(x,μ,Σ) = -length(μ)/2 * log(2π) - 1/2 * logdet(Σ) - 1/2 * (x-μ)' * inv(Σ) * (x-μ)

function fDrawInitialStatesBSPF(θ)

    Ψ0, Ψ1, Ψu, Φ1, Φε, vYlabs, vXlabs = fGetStateSpaceRep(θ)



    # Compute uconditional variance of states based on linearized model:
    nSlSS       = size(Φ1,1)
    P00vec      = myInv( diagm(vec(ones(nSlSS^2,1))) - kron(Φ1,Φ1) ) * vec(Φε * Φε')
    P00         = reshape(P00vec,(nSlSS,nSlSS))

    vslSS       = fMyMNDraw(zeros(nSlSS),P00)

    k_hat       = vslSS[vXlabs.=="k'"][1]
    lhat        = vslSS[vXlabs.=="l"][1]
    yhat        = vslSS[vXlabs.=="y"][1]
    ihat        = vslSS[vXlabs.=="i"][1]

    z           = vslSS[vXlabs.=="z"][1]
    b           = vslSS[vXlabs.=="b"][1]

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)

    K_          = kSS * (1+k_hat)
    L           = lSS * (1+lhat)
    Y           = ySS * (1+yhat)
    I           = iSS * (1+ihat)

    Z           = Zstar * exp(z)
    B           = Bstar * exp(b)

    vs          = [K_,L,Y,I]
    vε          = [z,b]


    return vs, vε

end

function fBSPF(mData,Mbspf,nS,nε,θ,fΨ,fΦ,fFε,fDrawInitialStates,fSol_K_)

    # Bootstrap Particle Filter for linear state space of the form
    #           s_t  =  fΦ(s_{t-1},ε_t;θ)         ,   ε_t ∼ fFε(ε_{t-1};θ)
    #           y_t  =  fΨ(s_t;θ)          +  u_t ,   u_t ∼ N(0,Σ_u(θ))
    # Note: with slight abuse of notation, ε are not the shocks/innovations to the exogenous shock processes, but the exogenous processes themselves (e.g. TFP process z that evolves acc. to AR(1))
    #
    # Inputs:
    # - (T x nY) matrix of data for likelihood evaluation,
    # - number of particles, number of endogenous states, number of exogenous states
    # - vector of model parameters
    # - SS objects (functions)
    # - function to draw initial states (vs,vε)
    # - policy function for capital (because fΦ uses it as an argument): it encodes model solution for this model
    #
    # Output: loglik, array of particles, matrix of weights
    #


    ### ------------------------ SET UP ----------------------------------------

    T,nY        = size(mData)

    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)

    Σ_u         = diagm([σy,σi,σl].^2)


    ### ------------------------ INITIALIZE OBJECTS ----------------------------

    aParticlesEndoStates    = zeros(T+1,Mbspf,nS)
    aParticlesExoStates     = zeros(T+1,Mbspf,nε)
    mWeights                = zeros(T+1,Mbspf)
    vLogLik                 = zeros(T)



    ### ------------------------ INITIALIZE ALGORITHM (RECURSIONS) -------------

    for jj = 1:Mbspf

        vs, vε      = fDrawInitialStatesBSPF(θ)

        aParticlesEndoStates[1,jj,:] = vs

        aParticlesExoStates[1,jj,:] = vε

    end

    mWeights[1,:]       .= 1


    # vVarHere = aParticlesEndoStates[1,:,4]
    # histogram(vVarHere)

    # vVarHere = aParticlesExoStates[1,:,2]
    # histogram(vVarHere)



    ### ------------------------ RECURSIONS ------------------------------------

    for t = 1:T


        # --- a1) Forecasting ε_t : -------------------------------------------------
        # draw from p(ε_t|ε_{t-1},θ) = fFε(ε_{t-1},θ)

        vvε                 = [ fFε(aParticlesExoStates[t,jj,:],ρz,σz,ρb,σb) for jj = 1:Mbspf ]

        # vVarHere = getindex.(vvε,2)
        # histogram(vVarHere)


        # --- a2) Forecasting s_t : -------------------------------------------------
        # compute s_t = fΦ(s_{t-1},ε_t;θ)

        vvs                     = [ fΦ(aParticlesEndoStates[t,jj,:],vvε[jj],fSol_K_,θ,Zstar,Bstar) for jj = 1:Mbspf ]

        # vVarHere = getindex.(vvs,4)
        # histogram(vVarHere)


        # --- b) Forecasting y_t : -------------------------------------------------
        # evaluate p(y_t|s_t,θ), y_t = fΨ(s_t;θ) + u_t ,   u_t ∼ N(0,Σ_u(θ))

        vLogIncrWeights         = [ ll( mData[t,:], fΨ(vvs[jj]) , Σ_u )  for jj = 1:Mbspf ]


        # Likelihood increment:
        vHelp                   = vLogIncrWeights .+ log.(mWeights[t,:]) #log of w-tilde * W
        c                       = maximum(vHelp)
        vHelp                   = vHelp .- c

        vLogLik[t]              = log(mean(exp.(vHelp))) + c


        # --- c) Normalized Weights : ----------------------------------------------

        vNormWeights            = exp.(vHelp) ./ mean(exp.(vHelp))


        # --- d) Selection : --------------------------------------------------------

        ESS                     = Mbspf / mean( vNormWeights.^2 )
        doSelect                = ESS < Mbspf/2

        if doSelect == 1

            #vSelectionIndices       = rand(Categorical(vNormWeights./Mbspf),Mbspf,1)
            vSelectionIndices       = resampleNYFED(vNormWeights)
            [aParticlesEndoStates[t+1,jj,:]   = vvs[vSelectionIndices[jj]] for jj = 1:Mbspf]
            [aParticlesExoStates[t+1,jj,:]   = vvε[vSelectionIndices[jj]] for jj = 1:Mbspf]
            mWeights[t+1,:]         .= 1

        elseif doSelect == 0

            mWeights[t+1,:]         = vNormWeights
            [aParticlesEndoStates[t+1,jj,:]   = vvs[jj] for jj = 1:Mbspf]
            [aParticlesExoStates[t+1,jj,:]   = vvε[jj] for jj = 1:Mbspf]

        end

    end


    ### ------------------------ COMPUTE LL APPROX. ----------------------------

    loglik          = sum( vLogLik )


    return loglik, aParticlesEndoStates, aParticlesExoStates, mWeights, vLogLik

end


# Function to simulate data based on this nonlinear SS rep:

function fSimulate_mVFI(fΨ,fΦ,fFε,θ,fSol_K_,nY,nS,nε, T, s0=0, ε0=0)

    # Simulates Linear State Space of the form:
    #           s_t  =  fΦ(s_{t-1},ε_t;θ)         ,   ε_t ∼ fFε(ε_{t-1};θ)
    #           y_t  =  fΨ(s_t;θ)          +  u_t ,   u_t ∼ N(0,Σ_u(θ))


    ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

    cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)

    Σ_u_sqrt       = diagm([σy,σi,σl])
    nη             = 3


    mε      = zeros(nε,T+1)
    mS      = zeros(nS,T+1)
    mY      = zeros(nY,T)

    s0 == 0 ? s0 = [kSS,lSS,ySS,iSS] : nothing
    ε0 == 0 ? ε0 = zeros(nε,1) : nothing

    mS[:,1] = s0

    mu      = Σ_u_sqrt * randn(nη,T)

    for tt = 1:T
        mε[:,tt+1] = fFε(mε[:,tt],ρz,σz,ρb,σb)
        mS[:,tt+1] = fΦ(mS[:,tt],mε[:,tt+1],fSol_K_,θ,Zstar,Bstar)
        mY[:,tt]   = fΨ(mS[:,tt+1]) + mu[:,tt]
    end


    return mY', mS', mε', mu'

end



# -- Resulting Likelihood Function ---------------------------------------------


function fLL_VFI(θ)

    fSol_K_ = fSolveModelVFI(θ)

    Mbspf   = 5000

    nS      = 4
    nε      = 2

    return fBSPF(mData,Mbspf,nS,nε,θ,fΨ,fΦ,fFε,fDrawInitialStatesBSPF,fSol_K_)[1]

end


# Only used for assessing BSPF accuracy (in "main_RBC_CheckBSPF.jl"):

function fLL_VFI(θ,Mbspf)

    timeSolution = @elapsed fSol_K_ = fSolveModelVFI(θ)

    nS      = 4
    nε      = 2

    timePF = @elapsed loglik = fBSPF(mData,Mbspf,nS,nε,θ,fΨ,fΦ,fFε,fDrawInitialStatesBSPF,fSol_K_)[1]

    return loglik, timeSolution, timePF

end
