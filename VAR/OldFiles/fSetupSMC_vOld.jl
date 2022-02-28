
# Marko Mlikota, mlikota@sas.upenn.edu
# November 2021

# -------------------------------------------------------------------------------


# Sets up functions for Sequential Monte Carlo (SMC) algorithm.


# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !


# -------------------------------------------------------------------------------



# Before running the function fSMC, need to define:

#       - @everywhere fLLtilde(θ)
#       - @everywhere fLL(θ)
#       - @everywhere fPriorDraw() (only if not all initial particles are supplied under MT)
#       - @everywhere fIsDrawValid(θ)
#       - @everywhere fPriorLogEval(θ)

#       - @everywhere vParLabs
#       - sOutputPath (should end with /)
#       - sOutputFilePrefix


#       For more details on software structuring for SMC estimation, see "HowToSetupSMC-EstimationSoftware.txt"

#       Note: it's important to write functions that evaluate prior, LL and LL-tilde directly as logs, not first the actual pdf value and only then logs. This avoids getting a zero, which then gives -Inf. This becomes relevant when proposal and target are very far from each other.


# fSMC will return several objects (see below, in function) and also write the following (if enabled):

#       - FinalStats.csv : 1 x 4 array with #Stages (Nϕ), #Particles (N), runtime, logMDD
#       - StageAll_Stats.csv : Nϕ x 4 array with tempering parameter, mutation step scaling constant, rejection rate, log MDD increment
#       - Stage1_Particles.csv : N x nP array with actual particles
#       - Stage1_WLLP.csv : N x 4 array with weight, LL, LLtilde, prior for each particle
#       - Stage1_ParticleStats.csv : nP x 4 array with mean, std, 5th and 95th percentile for each parameter
#       - StageLast_Particles.csv, StageLast_WLLP.csv, StageLast_ParticleStats.csv : analogous to above



# -------------------------------------------------------------------------------

# Improvment suggestions:

# - adapt the functions that pertain to blocking based on the more efficient code written for VAR-ZLB:
    # thetax                  = zeros(n_theta)
    #
    # mRandomDraws            = Random.rand(MvNormal(zeros(nVary),Matrix(1.0I, nVary, nVary)))
    #
    # thetax[vParamVaries]    = thetaMut[vParamVaries] + cSMC* mPropVar_sqrt* mRandomDraws
    #
    # thetax[vParamVaries.==0] = thetaMut[vParamVaries.==0]

# - have SMC return struct, not list of objects...

# -------------------------------------------------------------------------------


@everywhere begin

    function fLegitProposalDraw(θ0,c,Σ)

        # Draws proposal for which log-likelihood(s) and prior can be evaluated

        while true

            v       = rand( MultivariateNormal(θ0, c^2*Σ) )

            fIsDrawValid(v) == 1 ? nothing : continue

            llv, lltv, pv =
            try
                fLL(v), fLLtilde(v), fPriorLogEval(v)
            catch
                99, 99, 99
            end

            llv != 99 && lltv != 99 && pv != 99 && ~isnan(llv) && ~isnan(lltv)  && ~isnan(pv) ? (return v, [llv, lltv], pv) : nothing
            # Note: -Inf cases are handled by fRWMH1 (it's just a bad draw; no need to restrict)

        end

    end

    function fLegitProposalDraw(θ0,c,Σ,vInd)

        # Draws proposal for the parameters specified by vInd, s.t. the log-likelihood(s) and prior can be evaluated for the overall draw (θ0 with new proposal draw for specified parameters)

        while true

            vb        = rand( MultivariateNormal(θ0[vInd], c^2*Σ[vInd,vInd]) )
            v         = copy(θ0)
            #v[setdiff(1:vP,vInd)] = θ0[setdiff(1:vP,vInd)]
            v[vInd]   = vb

            fIsDrawValid(v) == 1 ? nothing : continue

            llv, lltv, pv =
            try
                fLL(v), fLLtilde(v), fPriorLogEval(v)
            catch
                99, 99, 99
            end

            llv != 99 && lltv != 99 && pv != 99 && ~isnan(llv) && ~isnan(lltv)  && ~isnan(pv) ? (return v, [llv, lltv], pv) : nothing
            # Note: -Inf cases are handled by fRWMH1 (it's just a bad draw; no need to restrict)

        end

    end

    # θ0,vllθ0,pθ0,c,Σ,ϕn = mDraws[ind-1,:],mLogLiks[ind-1,:],vPriors[ind-1],c,Σ,ϕn
    # Random.seed!((432 + 100*iiPhi + 432*ii)*setSeed)
    function fRWMH1(θ0,vllθ0,pθ0,c,Σ,ϕn)

        # Propagates particle θ0 through one iteration of RWMH
        # Input: θ0, 2-element vector of log-likelihoods at θ0 (set first to zero for likelihood tempering), prior at θ0, proposal distribution objects c & Σ, exponent ϕn for acceptance probability
        # Output:  resulting particle, Boolean for whether proposal was accepted, vector of log-likelihoods at particle, prior at particle

        v, vllv , pv        = fLegitProposalDraw(θ0,c,Σ)
        llθ0, lltθ0         = vllθ0
        llv, lltv           = vllv
        ratio               = ϕn*llv + (1-ϕn)*lltv + pv - ϕn*llθ0 - (1-ϕn)*lltθ0 - pθ0
        # isnan(ratio) ? ratio = ϕn*(llv - llθ0) + (1-ϕn)*(lltv-lltθ0) : nothing
        if ϕn==1
            ratio = ϕn*(llv - llθ0) + pv - pθ0
        end
        α                   = max(0,min(1,exp(ratio))) #min(1,exp(ratio))
        accept              = rand( Bernoulli(α) )

        accept == 1 ? (return v, accept, vllv, pv) : (return θ0, accept, vllθ0, pθ0)

    end

    function fRWMH1(θ0,vllθ0,pθ0,c,Σ,ϕn,nB,vOrder)

        # Propagates particle θ0 through one iteration of RWMH using nB blocks
        # Input: θ0, 2-element vector of log-likelihoods at θ0 (set first to zero for likelihood tempering), prior at θ0, proposal distribution objects c & Σ, exponent ϕn for acceptance probability, number of blocks, vector of indices (e.g. randomized) to set up blocks
        # Output:  resulting particle, Boolean for whether proposal was accepted, vector of log-likelihoods at particle, prior at particle

        blockLength     = Int(floor(length(θ0)/nB))
        vAccept         = zeros(nB,1)

        for b = 1:nB

            b < nB ?  vIndices = (b-1)*blockLength .+ (1:blockLength) : vIndices = (b-1)*blockLength .+ (1:length(vOrder)-(b-1)*blockLength)
            vInd                = vOrder[vIndices]

            v, vllv , pv        = fLegitProposalDraw(θ0,c,Σ,vInd)
            llθ0, lltθ0         = vllθ0
            llv, lltv           = vllv
            ratio               = ϕn*llv + (1-ϕn)*lltv + pv - ϕn*llθ0 - (1-ϕn)*lltθ0 - pθ0
            #isnan(ratio) ? ratio = ϕn*(llv - llθ0) + (1-ϕn)*(lltv-lltθ0) : nothing
            if ϕn==1
                ratio           = ϕn*(llv - llθ0) + pv - pθ0
            end
            # println(isfinite(ratio))
            α                   = max(0,min(1,exp(ratio))) #min(1,exp(ratio))
            vAccept[b]          = rand( Bernoulli(α) )

            if vAccept[b] == 1
                θ0         = v
                vllθ0      = vllv
                pθ0        = pv
            end

            b == nB ? (return θ0, mean(vAccept), vllθ0, pθ0) : nothing

        end

    end


    # θ0,vllθ0,pθ0,c,Σ,ϕn,nDraws = mProposals[ii,:],mLogLiks[ii,:],vPriors[ii],c,Σ,ϕn,Nmh
    # Random.seed!((432 + 100*iiPhi + 432*ii)*setSeed)
    function fRWMH(θ0,vllθ0,pθ0,c,Σ,ϕn,nDraws)

        # Propagates θ0 through nDraws iterations of RWMH. See fRWMH1().

        mDraws          = Array{Float64}(undef,nDraws+1,length(θ0))
        vAccept         = Array{Float64}(undef,nDraws+1,1)
        mLogLiks        = Array{Float64}(undef,nDraws+1,2)
        vPriors         = Array{Float64}(undef,nDraws+1,1)

        mDraws[1,:]     = θ0
        mLogLiks[1,:]   = vllθ0
        vPriors[1]      = pθ0

        for i = 1:nDraws

            ind     = i + 1

            mDraws[ind,:], vAccept[ind], mLogLiks[ind,:], vPriors[ind] = fRWMH1(mDraws[ind-1,:],mLogLiks[ind-1,:],vPriors[ind-1],c,Σ,ϕn)

        end

        return mDraws[end,:], sum(vAccept[2:end].==0)/nDraws, mLogLiks[end,:], vPriors[end]

    end

    function fRWMH(θ0,vllθ0,pθ0,c,Σ,ϕn,nDraws,nB,vOrder)

        # Propagates θ0 through nDraws iterations of RWMH. See fRWMH1().

        mDraws          = Array{Float64}(undef,nDraws+1,length(θ0))
        vAccept         = Array{Float64}(undef,nDraws+1,1)
        mLogLiks        = Array{Float64}(undef,nDraws+1,2)
        vPriors         = Array{Float64}(undef,nDraws+1,1)

        mDraws[1,:]     = θ0
        mLogLiks[1,:]   = vllθ0
        vPriors[1]      = pθ0

        for i = 1:nDraws

            ind     = i + 1

            mDraws[ind,:], vAccept[ind], mLogLiks[ind,:], vPriors[ind] = fRWMH1(mDraws[ind-1,:],mLogLiks[ind-1,:],vPriors[ind-1],c,Σ,ϕn,nB,vOrder)

        end

        return mDraws[end,:], 1-sum(vAccept[2:end])/nDraws, mLogLiks[end,:], vPriors[end]

    end

    function fLegitPriorDraw()

        # Generates valid prior draws

        for ii = 1:1e10

            v       = fPriorDraw()
            fIsDrawValid(v) == 1 ? nothing : continue

            llv, lltv, pv =
            try
                fLL(v), fLLtilde(v), fPriorLogEval(v)
            catch
                99, 99, 99
            end

            llv != 99 && lltv != 99 && pv != 99 && ~isnan(llv) && ~isnan(lltv)  && ~isnan(pv) &&
            llv != -Inf && lltv != -Inf && pv != -Inf ? (return v, [llv, lltv], pv, ii) : nothing
            # Note: here actually need to restrict llv,lltv,pv to be finite because we need valid initial draws (i.e. probability under prior or either likelihood should not be 0).

        end

    end

    function fLegitPriorDrawPartial(θ0,vPIndDraw)

        # Generates valid prior draws for parameters indexed in vPIndDraw (Boolean), taking rest from θ0

        v                  = zeros(length(θ0))

        v[vPIndDraw.==0]   = θ0[vPIndDraw.==0]


        for ii =1:1e10

            vAllDraws         = fPriorDraw()
            v[vPIndDraw]      = vAllDraws[vPIndDraw]
            fIsDrawValid(v) == 1 ? nothing : continue

            llv, lltv, pv =
            try
                fLL(v), fLLtilde(v), fPriorLogEval(v)
            catch
                99, 99, 99
            end

            llv != 99 && lltv != 99 && pv != 99 && ~isnan(llv) && ~isnan(lltv)  && ~isnan(pv) &&
            llv != -Inf && lltv != -Inf && pv != -Inf ? (return v, [llv, lltv], pv, ii) : nothing
            # Note: here actually need to restrict llv,lltv,pv to be finite because we need valid initial draws (i.e. probability under prior or either likelihood should not be 0).

        end

    end


end


function fESS(ϕ,ϕOld,mLL,vOldWeights)

    # Computes Effective Sample Size (ESS), vector of normalized weights and MDD, given
    # - new and old ϕ,
    # - Nx2-matrix of log-likelihoods ("true",then tilde in rows)
    # - vector of old weights.

    N                   = size(mLL)[1]

    exponent            = ϕ-ϕOld
    vDiffInLogLiks      = mLL[:,1]-mLL[:,2]
    maxLL               = maximum(vDiffInLogLiks)
    vIncrWeights        = exp.( exponent * ( vDiffInLogLiks .- maxLL ) )

    vNewWeights         = vIncrWeights.* vOldWeights

    meanNewWeights      = mean(vNewWeights)
    vNormWeights        = vNewWeights ./meanNewWeights

    ESS                 = N^2 / sum(vNormWeights.^2)

    logMDD              = log(meanNewWeights) + exponent*maxLL

    return ESS, vNormWeights, logMDD

end

function fFindBounds(ϕOld,findϕ,stepsize)

    while true

        findϕ(ϕOld) >= 0 ? ϕOld  += stepsize : return ϕOld

    end

end


function fGetParticleStats(vWeights,mParticles,vQuantiles=[0.05,0.95])

    nP              = size(mParticles,2)
    nQ              = length(vQuantiles)

    vProbWeights    = ProbabilityWeights(vWeights)
    vMeans          = mean(mParticles,vProbWeights,dims=1)
    vVars           = diag( cov(mParticles,vProbWeights) )
    vStds           = sqrt.(vVars)
    mCI             = zeros(nP,nQ)
    [mCI[pp,:]      = quantile(mParticles[:,pp],vProbWeights,vQuantiles) for pp = 1:nP]

    mPostStats      = [vMeans' vStds mCI]

    return mPostStats

end


function fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed=0,mInitParticles=zeros(10,5),vPIndDraw=zeros(5),writeFiles=1)

    # Inputs:
    # - tSettings = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,vP,c0,accStar,Nbar,showMessages), where vP is number of particles that vary (I use tuple to combine different types; int vs float)
    # - sOutputPath and sOutputFilePrefix are strings
    # - setSeed is integer (when zero, no seed is set)
    # - matrix with initial particles (e.g. for model tempering or generalized data tempering); dimensions are N x vP
    # - vector of length vP indicating whether the corresponding parameter needs to be drawn from the prior or can be taken from mInitParticles


    timeInit        = time_ns()

    ### ------------------------ SET UP --------------------------------------------

    N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,vP,c0,accStar,Nbar,showMessages = tSettings

    if useAdaptive == 0
        vϕ          = [ϕ_Nϕ* (n/Nϕ)^λ for n = 0:Nϕ]
    elseif useAdaptive == 1
        Nϕ          = 5000
        vϕ          = zeros(Nϕ+1,1)
    end

    fAdjMStep(x)    = 0.95 + 0.1* exp(16*(x-accStar))/(1 + exp(16*(x-accStar)) )


    ### ------------------------ INITIALIZE OBJECTS --------------------------------

    aParticles      = SharedArray{Float64,3}((N,vP,Nϕ+1))
    mWeights        = Array{Float64}(undef,N,Nϕ+1)
    aLogLiks        = SharedArray{Float64,3}((N,2,Nϕ+1))
    mPriors         = SharedArray{Float64,2}((N,Nϕ+1))
    mRejRate        = SharedArray{Float64,2}((N,Nϕ+1))
    vRejRate        = Array{Float64}(undef,Nϕ+1,1)
    vSstep          = Array{Int}(undef,Nϕ+1,1)
    vSstep[1]       = 0
    vESS            = Array{Float64}(undef,Nϕ+1,1)
    vc              = Array{Float64}(undef,Nϕ+1,1)
    vLogMDD         = Array{Float64}(undef,Nϕ+1,1)

    @eval @everywhere N         = $N
    @eval @everywhere Nmh       = $Nmh
    @eval @everywhere setSeed   = $setSeed


    ### ------------------------ INITIALIZE ALGORITHM (RECURSIONS) -------------


    mWeights[:,1]           .= 1
    vESS[1]                 = N^2 / sum(mWeights[:,1].^2)

    thereAreInitParticles   = sum(mInitParticles[1,:]) != 0
    needToDrawMore          = sum(vPIndDraw) > 0    # only relevant if above gives true


    if thereAreInitParticles

        aParticles[:,:,1]   = mInitParticles

        if needToDrawMore

            mParticlesStage1    = aParticles[:,:,1]
            @eval @everywhere mParticlesStage1   = $mParticlesStage1
            @eval @everywhere vPIndDraw          = $vPIndDraw

            mPriorParticles     = SharedArray{Float64,2}((N,vP))
            vNDrawsNeeded       = SharedArray{Float64,2}((N,1))

            @sync @distributed for ii = 1:N

                setSeed != 0 ? Random.seed!(100*ii*setSeed) : nothing # because LL evaluation might be stochastic
                θ               = mParticlesStage1[ii,:]
                mPriorParticles[ii,:], aLogLiks[ii,:,1], mPriors[ii,1], vNDrawsNeeded[ii] = fLegitPriorDrawPartial(θ,vPIndDraw)

            end

            aParticles[:,:,1]   = mPriorParticles

            probPriorDrawValid  = N/sum(vNDrawsNeeded)
            vLogMDD[1]          = -log(probPriorDrawValid)

        else

            mParticlesStage1    = aParticles[:,:,1]
            @eval @everywhere mParticlesStage1   = $mParticlesStage1

            @sync @distributed for ii = 1:N

                setSeed != 0 ? Random.seed!(100*ii*setSeed) : nothing # because LL evaluation might be stochastic
                θ                   = mParticlesStage1[ii,:]

                aLogLiks[ii,1,1], aLogLiks[ii,2,1], mPriors[ii,1] =
                try
                    fLL(θ), fLLtilde(θ), fPriorLogEval(θ)
                catch
                    99, 99, 99
                end

            end
            NDidntWork          = sum(aLogLiks[:,1,1].==99.0) + sum(aLogLiks[:,1,1].==-Inf) + sum(isnan.(aLogLiks[:,1,1]))

            NDidntWork > 0 ? throw(string("ERROR: Number of particles where initialization didn't work: ",NDidntWork)) : nothing

            vLogMDD[1]          = 0

        end

    else                                    ### DRAW FROM PRIOR

        mPriorParticles         = SharedArray{Float64,2}((N,vP))
        vNDrawsNeeded           = SharedArray{Float64,2}((N,1))

        @sync @distributed for ii = 1:N
            #println(ii)
            setSeed != 0 ? Random.seed!(100*ii*setSeed) : nothing
            if vP == 1
                mPriorParticles[ii,1], aLogLiks[ii,:,1], mPriors[ii,1], vNDrawsNeeded[ii] = fLegitPriorDraw()
            else
                mPriorParticles[ii,:], aLogLiks[ii,:,1], mPriors[ii,1], vNDrawsNeeded[ii] = fLegitPriorDraw()
            end
        end

        aParticles[:,:,1]       = mPriorParticles

        probPriorDrawValid      = N/sum(vNDrawsNeeded)
        vLogMDD[1]              = -log(probPriorDrawValid)

    end


    showMessages == 1 ? println("") : nothing
    showMessages == 1 ? println("Initialization Done") : nothing


    # Write stats for initial stage:

    if writeFiles == 1

        savename            =  sOutputFilePrefix * "Stage1_Particles.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(vParLabs)), header=false)
        CSV.write(sOutputPath * savename,  DataFrame(aParticles[:,:,1]), append=true)

        savename            =  sOutputFilePrefix * "Stage1_WLLP.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(["weight","ll","llt","prior"])), header=false)
        CSV.write(sOutputPath * savename,  DataFrame([mWeights[:,1] aLogLiks[:,:,1] mPriors[:,1]]), append=true)

        mParticleStats = fGetParticleStats(mWeights[:,1],aParticles[:,:,1])
        savename            =  sOutputFilePrefix * "Stage1_ParticleStats.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(["mean","std","lower","upper"])), header=false)
        CSV.write(sOutputPath * savename,  DataFrame(mParticleStats), append=true)

    end



    ### ------------------------ C-S-M ITERATIONS ----------------------------------

    for iiPhi = 1:Nϕ


        useAdaptive == 1 && vϕ[iiPhi] >= ϕ_Nϕ ? break : nothing

        showMessages == 1 ? println("") : nothing
        showMessages == 1 ? println("Iteration: n = ", iiPhi) : nothing

        nPhi = iiPhi+1
        @eval @everywhere nPhi = $nPhi


        # --- 1. Correction ---

        if useAdaptive == 0

            ESS, vNormWeights, logMDD = fESS(vϕ[nPhi],vϕ[nPhi-1],aLogLiks[:,:,nPhi-1],mWeights[:,nPhi-1])

        else

            findϕ(ϕ)       = fESS(ϕ,vϕ[nPhi-1],aLogLiks[:,:,nPhi-1],mWeights[:,nPhi-1])[1] - α*vESS[nPhi-1]

            if findϕ(ϕ_Nϕ) >= 0
                vϕ[nPhi]      = ϕ_Nϕ
            else
                stepsize   = 0.001
                ϕProp      = fFindBounds(vϕ[nPhi-1],findϕ,stepsize)

                vϕ[nPhi]   = find_zero(findϕ,(ϕProp-stepsize,ϕProp))
            end

            ESS, vNormWeights, logMDD = fESS(vϕ[nPhi],vϕ[nPhi-1],aLogLiks[:,:,nPhi-1],mWeights[:,nPhi-1])

        end

        vLogMDD[nPhi]      = logMDD #log( sum(vIncrWeights.*mWeights[:,nPhi-1])/N )

        showMessages == 1 ? println("- C done:  ϕ = ",round(vϕ[nPhi]; digits=4)) : nothing


        # --- 2. Selection ---

        ρn                  = ESS < Nbar || vϕ[nPhi] >= ϕ_Nϕ
        vSstep[nPhi]        = ρn
        vESS[nPhi]          = ESS

        if ρn == 1

            setSeed != 0 ? Random.seed!((432 + 100*iiPhi)*setSeed) : nothing

            vIndices            = resampleNYFED(vNormWeights) #rand(Categorical(vNormWeights./N),N,1)
            mProposals          = reshape( aParticles[vIndices,:,nPhi-1], N, :)
            mWeights[:,nPhi]    = vec(ones(N,1))
            mLogLiks            = reshape(aLogLiks[vIndices,:,nPhi-1],N,2)
            vPriors             = mPriors[vIndices,nPhi-1]
            vESS[nPhi]          = N

        elseif ρn == 0

            mProposals          = aParticles[:,:,nPhi-1]
            mWeights[:,nPhi]    = vNormWeights
            mLogLiks            = aLogLiks[:,:,nPhi-1]
            vPriors             = mPriors[:,nPhi-1]

        end

        showMessages == 1 && ρn == 1 ? println("- S done: ESS =",vESS[nPhi]) : nothing


        # --- 3. Mutation ---

        vw                      = ProbabilityWeights(vNormWeights)
        Σ                       = cov(aParticles[:,:,nPhi-1],vw)

        iiPhi == 1 ? c = c0 : c = fAdjMStep(1-vRejRate[nPhi-1])*vc[nPhi-1]
        vc[nPhi]                = c

        @eval @everywhere c     = $c
        @eval @everywhere Σ     = $Σ
        @eval @everywhere ϕn    = $vϕ[nPhi]

        if nB == 1

            mStepTime = @elapsed @sync @distributed for ii = 1:N

                setSeed != 0 ? Random.seed!((432 + 100*iiPhi + 432*ii)*setSeed) : nothing

                aParticles[ii,:,nPhi], mRejRate[ii,nPhi], aLogLiks[ii,:,nPhi], mPriors[ii,nPhi]  = fRWMH(mProposals[ii,:],mLogLiks[ii,:],vPriors[ii],c,Σ,ϕn,Nmh)

            end

        else

            mStepTime = @elapsed @sync @distributed for ii = 1:N

                setSeed != 0 ? Random.seed!((432 + 100*iiPhi + 432*ii)*setSeed) : nothing

                vOrder      = shuffle(1:vP)

                aParticles[ii,:,nPhi], mRejRate[ii,nPhi], aLogLiks[ii,:,nPhi], mPriors[ii,nPhi]  = fRWMH(mProposals[ii,:],mLogLiks[ii,:],vPriors[ii],c,Σ,ϕn,Nmh,nB,vOrder)

            end

        end

        vRejRate[nPhi]         = sum(mRejRate[:,nPhi])/N

        showMessages == 1 ? println("- M done:  c = ",round(c; digits=4),", AR = ",round(1-vRejRate[nPhi]; digits=4),", ",round(mStepTime;digits=2),"s") : nothing

    end


    useAdaptive == 1 ? Nϕ = findfirst(vϕ.>= ϕ_Nϕ)[1]-1 : nothing

    timeTotal      = signed(time_ns()-timeInit)/1000000000
    showMessages == 1 ? println("Total Time = ",timeTotal) : nothing


    # Write stats for last stage and some for all stages:

    if writeFiles == 1

        savename            =  sOutputFilePrefix * "StageLast_Particles.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(vParLabs)), header=false)
        CSV.write(sOutputPath * savename,  DataFrame(aParticles[:,:,Nϕ+1]), append=true)

        savename            =  sOutputFilePrefix * "StageLast_WLLP.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(["weight","ll","llt","prior"])), header=false)
        CSV.write(sOutputPath * savename,  DataFrame([mWeights[:,Nϕ+1] aLogLiks[:,:,Nϕ+1] mPriors[:,Nϕ+1]]), append=true)

        mParticleStats = fGetParticleStats(mWeights[:,Nϕ+1],aParticles[:,:,Nϕ+1])
        savename            =  sOutputFilePrefix * "StageLast_ParticleStats.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(["mean","std","lower","upper"])), header=false)
        CSV.write(sOutputPath * savename,  DataFrame(mParticleStats), append=true)

        savename            =  sOutputFilePrefix * "StageAll_Stats.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(["phi","c","rejRate","logMDD"])), header=false)
        CSV.write(sOutputPath * savename,  DataFrame([vϕ[1:Nϕ+1] vc[1:Nϕ+1] vRejRate[1:Nϕ+1] vLogMDD[1:Nϕ+1]]), append=true)

        savename            =  sOutputFilePrefix * "FinalStats.csv"
        CSV.write(sOutputPath * savename,  DataFrame(permutedims(["Nphi","N","time","logMDD"])), header=false)
        CSV.write(sOutputPath * savename,  DataFrame([Nϕ N timeTotal sum(vLogMDD[1:Nϕ+1])]), append=true)

    end


    return aParticles[:,:,1:Nϕ+1], mWeights[:,1:Nϕ+1], aLogLiks[:,1,1:Nϕ+1], mPriors[:,1:Nϕ+1], vRejRate[1:Nϕ+1], vSstep[1:Nϕ+1], vESS[1:Nϕ+1], vc[1:Nϕ+1], vϕ[1:Nϕ+1], vLogMDD[1:Nϕ+1], timeTotal

end
