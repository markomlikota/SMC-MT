# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Sets up functions for Sequential Monte Carlo (SMC) algorithm.


# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !



# -------------------------------------------------------------------------------


# Before running the function fSMC, need to define:
#
#       - @everywhere fLLtilde(θ)
#       - @everywhere fLL(θ)
#       - @everywhere fPriorDraw() (only if not all initial particles are supplied under MT)
#       - @everywhere fIsDrawValid(θ)
#       - @everywhere fPriorLogEval(θ)
#
#       - @everywhere vParLabs
#       - sOutputPath (should end with /)
#       - sOutputFilePrefix


# fSMC will return several objects (see below, in function) and also write the following (if enabled):
#
#       - FinalStats.csv :              1 x 4 array with #Stages (Nϕ), #Particles (N), runtime, logMDD
#       - StageAll_Stats.csv :          Nϕ x 4 array with tempering parameter, mutation step scaling constant, rejection rate, log MDD increment
#       - Stage1_Particles.csv :        N x nP array with actual particles
#       - Stage1_WLLP.csv :             N x 4 array with weight, LL, LLtilde, prior for each particle
#       - Stage1_ParticleStats.csv :    nP x 4 array with mean, std, 5th and 95th percentile for each parameter
#       - StageLast_Particles.csv, StageLast_WLLP.csv, StageLast_ParticleStats.csv : analogous to above


# Notes:
#
#   - it's important to write functions that evaluate prior, LL and LL-tilde directly as logs, not first the actual pdf value and only then logs.
#     This avoids getting a zero, which then gives -Inf. This becomes relevant when proposal and target are very far from each other.
#
#   - prior draw validity (as it affects variable probPriorDrawValid below) assessed not only via fIsDrawValid(θ),
#     but also via fLL(θ), fLLtilde(θ) and fPriorLogEval(θ) bcs we require them to be computable (≂̸99,≂̸NaN) and finite
#     for a prior draw to be valid (see fLegitPriorDraw() or fLegitPriorDrawPartial())
#     As opposed to that, proposal draw validity in mutation step is only assessed via fIsDrawValid(θ) and computability (≂̸99,≂̸NaN)
#     of the latter three functions (see fLegitProposalDraw()).
#     The case where one of these three functions gives -Inf is unlikely (impossible) though if one computes directly the log pdf (prior or LL)
#     instead of computing the pdf and then taking logs (in case pdf value is low, log could give -Inf)
#
#   - For more details on software structuring for SMC estimation, see "HowToSetupSMC-EstimationSoftware.txt"


# -------------------------------------------------------------------------------


@everywhere begin

    function fLegitProposalDraw(θ0,c,Σ,vDoesParamVary)

        # Draws proposal for which log-likelihood(s) and prior can be evaluated

        while true

            v                           = zeros(length(θ0))
            v[vDoesParamVary.==false]   = θ0[vDoesParamVary.==false]
            v[vDoesParamVary]           = rand( MultivariateNormal(θ0[vDoesParamVary], c^2*Σ) )

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


    function fRWMH1(θ0,vllθ0,pθ0,c,Σ,vDoesParamVary,ϕn)

        # Propagates particle θ0 through one iteration of RWMH
        # Input: θ0, 2-element vector of log-likelihoods at θ0 (set first to zero for likelihood tempering), prior at θ0, proposal distribution objects c & Σ, exponent ϕn for acceptance probability
        # Output:  resulting particle, Boolean for whether proposal was accepted, vector of log-likelihoods at particle, prior at particle

        v, vllv , pv        = fLegitProposalDraw(θ0,c,Σ,vDoesParamVary)
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

    function fRWMH(θ0,vllθ0,pθ0,c,Σ,vDoesParamVary,ϕn,nDraws)

        # Propagates θ0 through nDraws iterations of RWMH. See fRWMH1().

        mDraws          = Array{Float64}(undef,nDraws+1,length(θ0))
        vAccept         = Array{Float64}(undef,nDraws+1,1)
        mLogLiks        = Array{Float64}(undef,nDraws+1,2)
        vPriors         = Array{Float64}(undef,nDraws+1,1)

        mDraws[1,:]     = θ0
        mLogLiks[1,:]   = vllθ0
        vPriors[1]      = pθ0

        for iii = 1:nDraws

            ind     = iii + 1

            mDraws[ind,:], vAccept[ind], mLogLiks[ind,:], vPriors[ind] = fRWMH1(mDraws[ind-1,:],mLogLiks[ind-1,:],vPriors[ind-1],c,Σ,vDoesParamVary,ϕn)

        end

        return mDraws[end,:], sum(vAccept[2:end])/nDraws, mLogLiks[end,:], vPriors[end]

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


function fSMC(tSettings,sOutputPath,sOutputFilePrefix="",writeFiles=1,setSeed=0,vaMT=[],vaContinue=[])


    # Inputs:
    # - tSettings = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages), where nP is number of particles that vary (I use tuple to combine different types; int vs float)
    # - sOutputPath and sOutputFilePrefix are strings
    # - writeFiles: 0: no, 1: only first and last stage (and allstage-files), 2: every stage
    # - setSeed is integer (when zero, no seed is set)
    # - vaMT is 2x1 and contains the following:
    #   - matrix with initial particles (e.g. for model tempering or generalized data tempering); dimensions are N x nP
    #   - vector of length nP indicating whether the corresponding parameter needs to be drawn from the prior or can be taken from mInitParticles
    # - vaContinue is 3x1 and contains the following matrices:
    #   - mContinue_StageAll_Stats
    #   - mContinue_WLLP
    #   - mContinue_Particles
    #   These matrices are output from an interrupted SMC run. When supplied, SMC continues from where it left off.


    timeInit        = time_ns()



    ### ------------------------ SET UP --------------------------------------------


    N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages = tSettings

    if useAdaptive == 0
        vϕ          = [ϕ_Nϕ* (n/Nϕ)^λ for n = 0:Nϕ]
    elseif useAdaptive == 1
        vϕ           = Vector{Float64}(undef,1) # (yet) undefined length
        vϕ[1]        = 0
    end

    fAdjMStep(x)    = 0.95 + 0.1* exp(16*(x-accStar))/(1 + exp(16*(x-accStar)) )


    # Extract any objects inside vaMT and vaContinue:

    if length(vaMT) == 0 # i.e. when vaMT = [], i.e. not supplied
        mInitParticles = zeros(10,5)
        vPIndDraw      = zeros(5)
    else
        mInitParticles = vaMT[1]
        vPIndDraw      = vaMT[2]
    end

    if length(vaContinue) == 0 # i.e. when vaContinue = [], i.e. not supplied
        mContinue_StageAll_Stats    = 0
        mContinue_WLLP              = 0
        mContinue_Particles         = 0
    else
        mContinue_StageAll_Stats    = vaContinue[1]
        mContinue_WLLP              = vaContinue[2]
        mContinue_Particles         = vaContinue[3]
    end



    ### ------------------------ INITIALIZE OBJECTS --------------------------------


    vmParticles       = Vector{Matrix{Float64}}(undef,1)
    vmLogLiks         = Vector{Matrix{Float64}}(undef,1)
    vvPriors          = Vector{Vector{Float64}}(undef,1)
    vvWeights         = Vector{Vector{Float64}}(undef,1)
    vAcceptRate       = Vector{Float64}(undef,1)
    vSstep            = Vector{Float64}(undef,1)
    vESS              = Vector{Float64}(undef,1)
    vc                = Vector{Float64}(undef,1)
    vLogMDD           = Vector{Float64}(undef,1)
    vTimes            = Vector{Float64}(undef,1)
    # all set up this way because #Stages Nϕ is (yet) undefined under adaptive tempering

    @eval @everywhere N         = $N
    @eval @everywhere Nmh       = $Nmh
    @eval @everywhere setSeed   = $setSeed



    ### ------------------------ INITIALIZE ALGORITHM (RECURSIONS) -------------


    continueInterruptedRun  =  size(mContinue_StageAll_Stats,1) > 1

    if continueInterruptedRun

        nPhiInterrupted     = size(mContinue_StageAll_Stats,1)
        iiPhiInterrupted    = nPhiInterrupted -1

        # Some objects need to be re-created in full:
        vLogMDD[1]          = mContinue_StageAll_Stats[1,4]
        vTimes[1]           = mContinue_StageAll_Stats[1,5]
        for iiPhi = 1:iiPhiInterrupted
            nPhi = iiPhi+1
            push!(vLogMDD,mContinue_StageAll_Stats[nPhi,4])
            push!(vTimes,mContinue_StageAll_Stats[nPhi,5])
        end

        #Some only for last stage:

        vmParticles[1]      = zeros(N,nP)
        vmLogLiks[1]        = zeros(N,2)
        vvPriors[1]         = zeros(N)
        vvWeights[1]        = ones(N)
        vAcceptRate[1]      = accStar
        vSstep[1]           = 0
        vESS[1]             = N
        vc[1]               = c0

        for iiPhi = 1:(iiPhiInterrupted-1)
            nPhi = iiPhi + 1

            push!(vmParticles,zeros(N,nP))
            push!(vmLogLiks,zeros(N,2))
            push!(vvPriors,zeros(N))
            push!(vvWeights,ones(N))
            push!(vAcceptRate,accStar)
            push!(vSstep,0)
            push!(vESS,N)
            push!(vc,c0)

            if useAdaptive == 1
                push!(vϕ,0)
            end

        end

        push!(vmParticles,mContinue_Particles)
        push!(vmLogLiks,mContinue_WLLP[:,2:3])
        push!(vvPriors,mContinue_WLLP[:,4])
        push!(vvWeights,mContinue_WLLP[:,1])
        push!(vAcceptRate,mContinue_StageAll_Stats[end,3])
        push!(vSstep,0) #no need to recreate at all
        push!(vESS,mContinue_StageAll_Stats[end,6])
        push!(vc,mContinue_StageAll_Stats[end,2])
        if useAdaptive == 1
            push!(vϕ,mContinue_StageAll_Stats[end,1])
        end


        # Determine which particles are fixed:

        vDoesParamVary       = [length(unique(vmParticles[end][:,pp]))!=1 for pp = 1:nP]

        @eval @everywhere vDoesParamVary         = $vDoesParamVary


        # Add time needed to re-initialize things:

        timeThisStage       = signed(time_ns()-timeInit)/1000000000
        vTimes[end]         += timeThisStage


    else


        vvWeights[1]            = ones(N)
        vESS[1]                 = N^2 / sum(vvWeights[1].^2)
        vSstep[1]               = 0
        vAcceptRate[1]          = NaN

        thereAreInitParticles   = sum(mInitParticles[1,:]) != 0
        needToDrawMoreParams    = sum(vPIndDraw) > 0    # only relevant if above gives true


        if thereAreInitParticles

            if needToDrawMoreParams

                # Draw additional parameters, and combine with existing, external ones:

                @eval @everywhere mInitParticles   = $mInitParticles
                @eval @everywhere vPIndDraw        = $vPIndDraw

                mParticlesStage1    = SharedArray{Float64,2}((N,nP))
                vNDrawsNeeded       = SharedArray{Float64,1}((N))
                vPriorsStage1       = SharedArray{Float64,1}((N))
                mLogLiksStage1      = SharedArray{Float64,2}((N,2))

                @sync @distributed for ii = 1:N

                    setSeed != 0 ? Random.seed!(100*ii*setSeed) : nothing # because LL evaluation might be stochastic
                    θ               = mInitParticles[ii,:]
                    mParticlesStage1[ii,:], mLogLiksStage1[ii,:], vPriorsStage1[ii], vNDrawsNeeded[ii] = fLegitPriorDrawPartial(θ,vPIndDraw)

                end

                # Store particles that consist of combined parameter vectors:

                vmParticles[1]      = mParticlesStage1
                vmLogLiks[1]        = mLogLiksStage1
                vvPriors[1]         = vPriorsStage1

                probPriorDrawValid  = N/sum(vNDrawsNeeded)
                vLogMDD[1]          = -log(probPriorDrawValid)

            else

                # Verify that LL can be computed for existing particles:

                @eval @everywhere mInitParticles   = $mInitParticles

                vPriorsStage1       = SharedArray{Float64}((N))
                mLogLiksStage1      = SharedArray{Float64,2}((N,2))

                @sync @distributed for ii = 1:N

                    setSeed != 0 ? Random.seed!(100*ii*setSeed) : nothing # because LL evaluation might be stochastic
                    θ                   = mInitParticles[ii,:]

                    mLogLiksStage1[ii,1], mLogLiksStage1[ii,2], vPriorsStage1[ii] =
                    try
                        fLL(θ), fLLtilde(θ), fPriorLogEval(θ)
                    catch
                        99, 99, 99
                    end

                end

                NDidntWork          = sum(mLogLiksStage1[:,1].==99.0) + sum(mLogLiksStage1[:,1].==-Inf) + sum(isnan.(mLogLiksStage1[:,1]))

                NDidntWork > 0 ? throw(string("ERROR: Number of particles where initialization didn't work: ",NDidntWork)) : nothing


                # Store particles:

                vmParticles[1]      = mInitParticles
                vmLogLiks[1]        = mLogLiksStage1
                vvPriors[1]         = vPriorsStage1

                vLogMDD[1]          = 0

            end

        else

            # Draw particles from prior:

            mParticlesStage1        = SharedArray{Float64,2}((N,nP))
            vNDrawsNeeded           = SharedArray{Float64,1}((N))
            vPriorsStage1           = SharedArray{Float64,1}((N))
            mLogLiksStage1          = SharedArray{Float64,2}((N,2))

            @sync @distributed for ii = 1:N
                #println(ii)
                setSeed != 0 ? Random.seed!(100*ii*setSeed) : nothing
                if nP == 1
                    mParticlesStage1[ii,1], mLogLiksStage1[ii,:], vPriorsStage1[ii], vNDrawsNeeded[ii] = fLegitPriorDraw()
                else
                    mParticlesStage1[ii,:], mLogLiksStage1[ii,:], vPriorsStage1[ii], vNDrawsNeeded[ii] = fLegitPriorDraw()
                end
            end


            # Store particles:

            vmParticles[1]          = mParticlesStage1
            vmLogLiks[1]            = mLogLiksStage1
            vvPriors[1]             = vPriorsStage1

            probPriorDrawValid      = N/sum(vNDrawsNeeded)
            vLogMDD[1]              = -log(probPriorDrawValid)

        end

        showMessages == 1 ? println("") : nothing
        showMessages == 1 ? println("Initialization Done") : nothing


        # Write stats for initial stage:

        timeThisStage       = signed(time_ns()-timeInit)/1000000000
        vTimes[1]           = timeThisStage

        if writeFiles > 0

            savename            =  sOutputFilePrefix * "Stage1_Particles.csv"
            fMyCSVWRITE(sOutputPath * savename,vmParticles[1],vParLabs)

            savename            =  sOutputFilePrefix * "Stage1_WLLP.csv"
            fMyCSVWRITE(sOutputPath * savename,[vvWeights[1] vmLogLiks[1] vvPriors[1]],["weight","ll","llt","prior"])

            mParticleStats      = fGetParticleStats(vvWeights[1],vmParticles[1])
            savename            =  sOutputFilePrefix * "Stage1_ParticleStats.csv"
            fMyCSVWRITE(sOutputPath * savename,mParticleStats,["mean","std","lower","upper"])

            if writeFiles == 2

                try
                    mkdir(sOutputPath * "IntermediateStages")
                catch
                end

                savename            =  sOutputFilePrefix * "StageAll_Stats.csv"
                fMyCSVWRITE(sOutputPath * savename,[vϕ[1] vc[1] vAcceptRate[1] vLogMDD[1] vTimes[1] vESS[1]],["phi","c","accRate","logMDD","time","ESS"])

            end

        end


        # Determine which particles are fixed:

        vDoesParamVary       = [length(unique(vmParticles[1][:,pp]))!=1 for pp = 1:nP]

        @eval @everywhere vDoesParamVary         = $vDoesParamVary


    end



    ### ------------------------ C-S-M ITERATIONS ----------------------------------


    useAdaptive == 1 ? Nϕ = 1000000 : nothing

    continueInterruptedRun ? iiPhiFirst = (iiPhiInterrupted+1) : iiPhiFirst = 1

    for iiPhi = iiPhiFirst : Nϕ

        timeStartStage       = time_ns()

        useAdaptive == 1 && vϕ[iiPhi] >= ϕ_Nϕ ? break : nothing

        showMessages == 1 ? println("") : nothing
        showMessages == 1 ? println("Stage:     n = ", iiPhi) : nothing


        nPhi        = iiPhi+1           # iiPhi is actual stage i, goes from 0 to Nϕ (see analytical outline of algorithm),
        @eval @everywhere nPhi = $nPhi  # nPhi is helper to store objects (goes from 1 to Nϕ+1)


        # --- 1. Correction ---

        if useAdaptive == 0

            ESS, vNormWeights, logMDD = fESS(vϕ[nPhi],vϕ[nPhi-1],vmLogLiks[nPhi-1],vvWeights[nPhi-1])

        else

            findϕ(ϕ)       = fESS(ϕ,vϕ[nPhi-1],vmLogLiks[nPhi-1],vvWeights[nPhi-1])[1] - α*vESS[nPhi-1]

            if findϕ(ϕ_Nϕ) >= 0
                push!(vϕ,ϕ_Nϕ)
            else
                stepsize   = 0.001
                ϕProp      = fFindBounds(vϕ[nPhi-1],findϕ,stepsize)

                nextϕ      = find_zero(findϕ,(ϕProp-stepsize,ϕProp))
                push!(vϕ,nextϕ)
            end

            ESS, vNormWeights, logMDD = fESS(vϕ[nPhi],vϕ[nPhi-1],vmLogLiks[nPhi-1],vvWeights[nPhi-1])

        end

        push!(vLogMDD,logMDD)

        showMessages == 1 ? println("- C done:  ϕ = ",round(vϕ[nPhi]; digits=4)) : nothing


        # --- 2. Selection ---

        ρn                  = ESS < Nbar || vϕ[nPhi] >= ϕ_Nϕ
        push!(vSstep,ρn)

        if ρn == 1

            setSeed != 0 ? Random.seed!((432 + 100*iiPhi)*setSeed) : nothing

            vIndices            = resampleNYFED(vNormWeights) #rand(Categorical(vNormWeights./N),N,1)

            mProposalParticles          = reshape( vmParticles[nPhi-1][vIndices,:], N, : )
            mProposalLogLiks            = reshape( vmLogLiks[nPhi-1][vIndices,:], N, 2 )
            vProposalPriors             = vvPriors[nPhi-1][vIndices]

            push!(vESS,N)
            push!(vvWeights,ones(N))

        elseif ρn == 0

            mProposalParticles          = vmParticles[nPhi-1]
            mProposalLogLiks            = vmLogLiks[nPhi-1]
            vProposalPriors             = vvPriors[nPhi-1]

            push!(vESS,ESS)
            push!(vvWeights,vNormWeights)

        end

        showMessages == 1 && ρn == 1 ? println("- S done: ESS =",vESS[nPhi]) : nothing


        # --- 3. Mutation ---

        vw                      = ProbabilityWeights(vNormWeights)
        Σ                       = cov(vmParticles[nPhi-1][:,vDoesParamVary],vw)

        iiPhi == 1 ? c = c0 : c = fAdjMStep(vAcceptRate[nPhi-1]) * vc[nPhi-1]
        push!(vc,c)

        @eval @everywhere c     = $c
        @eval @everywhere Σ     = $Σ
        @eval @everywhere ϕn    = $vϕ[nPhi]

        mParticlesThisStage     = SharedArray{Float64,2}((N,nP))
        vPriorsThisStage        = SharedArray{Float64,1}((N))
        mLogLiksThisStage       = SharedArray{Float64,2}((N,2))
        vAcceptRateThisStage    = SharedArray{Float64,1}((N))

        if nB == 1

            mStepTime = @elapsed @sync @distributed for ii = 1:N

                setSeed != 0 ? Random.seed!((432 + 100*iiPhi + 432*ii)*setSeed) : nothing

                mParticlesThisStage[ii,:], vAcceptRateThisStage[ii], mLogLiksThisStage[ii,:], vPriorsThisStage[ii]  = fRWMH(mProposalParticles[ii,:],mProposalLogLiks[ii,:],vProposalPriors[ii],c,Σ,vDoesParamVary,ϕn,Nmh)

            end

        else

            throw("Blocking in M-step doesn't work yet")

        end

        push!(vmParticles,mParticlesThisStage)
        push!(vvPriors,vPriorsThisStage)
        push!(vmLogLiks,mLogLiksThisStage)

        push!(vAcceptRate,sum(vAcceptRateThisStage)/N)

        showMessages == 1 ? println("- M done:  c = ",round(c; digits=4),", AR = ",round(vAcceptRate[nPhi]; digits=4),", ",round(mStepTime;digits=2),"s") : nothing


        # --- 4. Write Files From This Iteration ---

        timeThisStage       = signed(time_ns()-timeStartStage)/1000000000
        # timeThisStage       = signed(time_ns()-timeInit)/1000000000 - sum(vTimes[1:(nPhi-1)])
        push!(vTimes,timeThisStage)

        if writeFiles == 2

            savename            =  sOutputFilePrefix * string("Stage",nPhi,"_Particles.csv")
            fMyCSVWRITE(sOutputPath * "IntermediateStages/" * savename,vmParticles[nPhi],vParLabs)

            savename            =  sOutputFilePrefix * string("Stage",nPhi,"_WLLP.csv")
            fMyCSVWRITE(sOutputPath * "IntermediateStages/" * savename,[vvWeights[end] vmLogLiks[end] vvPriors[end]],["weight","ll","llt","prior"])

            savename            =  sOutputFilePrefix * "StageAll_Stats.csv"
            fMyCSVAPPEND(sOutputPath * savename,[vϕ[nPhi] vc[nPhi] vAcceptRate[nPhi] vLogMDD[nPhi] vTimes[nPhi] vESS[nPhi]])

        end



    end


    useAdaptive == 1 ? Nϕ = findfirst(vϕ.>= ϕ_Nϕ)[1]-1 : nothing

    timeTotal      = sum(vTimes) #signed(time_ns()-timeInit)/1000000000
    showMessages == 1 ? println("Total Time = ",timeTotal) : nothing


    # Write stats for last stage and some for all stages:

    if writeFiles > 0

        savename            =  sOutputFilePrefix * "StageLast_Particles.csv"
        fMyCSVWRITE(sOutputPath * savename,vmParticles[end],vParLabs)

        savename            =  sOutputFilePrefix * "StageLast_WLLP.csv"
        fMyCSVWRITE(sOutputPath * savename,[vvWeights[end] vmLogLiks[end] vvPriors[end]],["weight","ll","llt","prior"])

        mParticleStats      = fGetParticleStats(vvWeights[end],vmParticles[end])
        savename            =  sOutputFilePrefix * "StageLast_ParticleStats.csv"
        fMyCSVWRITE(sOutputPath * savename,mParticleStats,["mean","std","lower","upper"])

        if writeFiles != 2
            savename            =  sOutputFilePrefix * "StageAll_Stats.csv"
            fMyCSVWRITE(sOutputPath * savename,[vϕ vc vAcceptRate vLogMDD vTimes vESS],["phi","c","accRate","logMDD","time","ESS"])
        end

        savename            =  sOutputFilePrefix * "FinalStats.csv"
        fMyCSVWRITE(sOutputPath * savename,[Nϕ N timeTotal sum(vLogMDD)],["Nphi","N","time","logMDD"])

    end


    return vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal

end
