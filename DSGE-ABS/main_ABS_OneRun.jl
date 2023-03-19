# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# August 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Main script to estimate RBC model using LT/MT.


# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec         = "1"     # data from which DGP to take?
Model1spec      = "L2"    # which model (LL) to estimate?
Prior1spec      = "L2"    # which prior to take?
SMC1spec        = "1"    # which values for SMC tuning parameters to take?

sTemperingType  = "LT"          # "LT" or "MT"; how to estimate M1?
ϕ_Nϕ            = 1.0           # stop early?


# Only relevant if sTemperingType == "MT":

Model0spec      = "L1"            # which model is M0?
Prior0spec      = "L1"           # which prior was taken for M0?
ϕ_Nϕ_M0         = 1.0             # early stopping for M0?
SMC0spec        = "1"    # which values for SMC tuning parameters were taken for M0?


continueInterruptedRun = 0 #if possible, continue from interrupted run (code determines whether last run was complete or not)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


nProcs = 8

using Distributed
addprocs(nProcs-length(procs()))


using Plots
pyplot()
using StatsBase
using DataFrames
using CSV
using DelimitedFiles
using StatsFuns

@everywhere using ParallelDataTransfer
@everywhere using PositiveFactorizations
@everywhere using Random
@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using LinearAlgebra
@everywhere using Interpolations
@everywhere using Roots
@everywhere using SolveDSGE



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


# @everywhere cd()
# @everywhere sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/RBC/")
@everywhere sMyPath            = pwd() * "/"


@everywhere include(sMyPath * "/Functions/fHelpersSMC.jl")
include(sMyPath * "/Functions/fSetupSMC.jl")
@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")
include(sMyPath * "/Functions/fFolderFileManagement.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------



# -- Set up directory for estimation output ------------------------------------


sSpecFolder             = string("dgp",DGPspec,"_m",Model1spec,"_pr",Prior1spec,"_mc",SMC1spec,"/")
sHelperFolder           = ""
if sTemperingType == "LT"
    sEstSpecFolder          = string( "LTphiLast" , Int(ϕ_Nϕ * 100) , "/" )
elseif sTemperingType == "MT"
    sEstSpecFolder          = string( "MTphiLast" , Int(ϕ_Nϕ_M0 * 100) , "_m", Model0spec, "_pr", Prior0spec, "_mc", SMC0spec, "/" )
end


sSpecPath       = sMyPath * "Output/" * sSpecFolder
sHelperPath     = sSpecPath * sHelperFolder
sEstSpecPath    = sHelperPath * sEstSpecFolder

try
    mkdir(sSpecPath)
catch
end

try
    mkdir(sHelperPath)
catch
end

try
    mkdir(sEstSpecPath)
catch
end



# -- Load Data, Model, Prior, SMC-info -----------------------------------------


# Load data:

sFileToRead         = "Data_" * string("dgp",DGPspec) * ".csv"
sPathToRead         = sMyPath * "Data/" * sFileToRead

# mData               = CSV.read(sPathToRead;header=1) # DOESN'T WORK YET ...
mData               = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)
@everywhere @eval mData     = $mData


# Load Model (Likelihood) and prior:

# This is now done inside "Priorspec"-functions, but keep it here just in case...:
# @everywhere @eval DGPspec        = $DGPspec # we need DGPspec, because prior fixes some parameters at true values, which are contained in DGPspec
# @everywhere include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )

@everywhere @eval DGPspec        = $DGPspec # we need DGPspec, because prior fixes some parameters at true values, which are contained in DGPspec


@everywhere @eval Model1spec     = $Model1spec
@everywhere @eval Prior1spec     = $Prior1spec

@everywhere include(sMyPath * "SpecFiles/" * string("script_Modelspec",Model1spec,".jl") )
@everywhere include(sMyPath * "SpecFiles/" * string("script_Priorspec",Prior1spec,".jl") )

if sTemperingType == "MT"
    @everywhere @eval Model0spec     = $Model0spec
    @everywhere @eval Prior0spec     = $Prior0spec

    @everywhere include(sMyPath * "SpecFiles/" * string("script_Modelspec",Model0spec,".jl") )
    @everywhere include(sMyPath * "SpecFiles/" * string("script_Priorspec",Prior0spec,".jl") )
end


# Load SMC tuning parameters:

include(sMyPath * "SpecFiles/" * string("script_SMCspec",SMC1spec,".jl") )


# -- Estimation ----------------------------------------------------------------


if sTemperingType == "LT"

    # Define LL, LL-tilde and prior functions :

    @everywhere begin

        #Need to use anonymous functions because they are changed every time Model or prior spec files are changed.

        function fLL(θ)
            fLLhelp     = eval(Meta.parse(string("fLL_",Model1spec)))
            return fLLhelp(θ)
        end
        fLLtilde(θ)     = 0
        function fPriorDraw()
            fPDhelp     = eval(Meta.parse(string("fPriorDraw_",Prior1spec)))
            return fPDhelp()
        end
        function fIsDrawValid(θ)
            fIDVhelp         =  eval(Meta.parse(string("fIsDrawValid_",Prior1spec)))
            return minimum(fIDVhelp(θ))
        end
        function fPriorLogEval(θ)
            fPLEhelp         =  eval(Meta.parse(string("fPriorLogEval_",Prior1spec)))
            return sum(fPLEhelp(θ))
        end

    end


    # Setup inputs to SMC algorithm:

    nP              = length(vParLabs)
    showMessages    = 1
    tSettings       = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages)
    @everywhere @eval tSettings     = $tSettings

    setSeed             = 547
    sOutputFilePrefix   = ""
    sOutputPath         = sEstSpecPath


    # Continue interrupted run?

    doesPreviousRunExist = fDoesFileExist(sOutputPath, sOutputFilePrefix* string("StageAll_Stats.csv"))

    if doesPreviousRunExist && continueInterruptedRun == 1

        mContinue_StageAll_Stats    = fGetFile(sOutputPath, sOutputFilePrefix* string("StageAll_Stats.csv"))

        previousRunIncomplete = mContinue_StageAll_Stats[end,1] < 1.0

        if previousRunIncomplete

            nPhiInterrupted             = size(mContinue_StageAll_Stats,1)

            mContinue_WLLP              = fGetFile(sOutputPath * "IntermediateStages/", sOutputFilePrefix* string("Stage",nPhiInterrupted,"_WLLP.csv"))

            mContinue_Particles         = fGetFile(sOutputPath * "IntermediateStages/", sOutputFilePrefix* string("Stage",nPhiInterrupted,"_Particles.csv"))

            # vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,2,setSeed,zeros(10,5),zeros(5),mContinue_StageAll_Stats,mContinue_WLLP,mContinue_Particles)

            vaMT       = []
            vaContinue = [mContinue_StageAll_Stats, mContinue_WLLP, mContinue_Particles]

            vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,2,setSeed,vaMT,vaContinue)

        end

    end

    if continueInterruptedRun == 0 || doesPreviousRunExist == false

        # Run SMC:

        vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,2,setSeed)

    end


elseif sTemperingType == "MT"

    # Define LL, LL-tilde and prior functions :

    # Note: compared to VAR-SV application, here we take vector of parameters used in SMC to be same across both (all) models and allow some of them to be fixed.
    # Also, things are not written in most general way; one could
    #   - have parameters that are fixed under M1 but vary under M0 (which one needs to fix then in prior under M1 at posterior mean under M0 (e.g.))
    #   - different prior distributions for parameters that vary under both models
    # We do, however, have a parameter that varies under "mVFI" but not under "mL1": ϕ2


    @everywhere begin

        # First "translate" model functions to M1 and M0:

        function fLL_M1(θ)
            fLLM1help     = eval(Meta.parse(string("fLL_",Model1spec)))
            return fLLM1help(θ)
        end
        function fPriorDraw_M1()
            fPDM1help     = eval(Meta.parse(string("fPriorDraw_",Prior1spec)))
            return fPDM1help()
        end
        function fIsDrawValid_M1(θ)
            fIDVM1help         =  eval(Meta.parse(string("fIsDrawValid_",Prior1spec)))
            return fIDVM1help(θ)
        end
        function fPriorLogEval_M1(θ)
            fPLEM1help         =  eval(Meta.parse(string("fPriorLogEval_",Prior1spec)))
            return fPLEM1help(θ)
        end


        function fLL_M0(θ)
            fLLM0help     = eval(Meta.parse(string("fLL_",Model0spec)))
            return fLLM0help(θ)
        end
        function fPriorDraw_M0()
            fPDM0help     = eval(Meta.parse(string("fPriorDraw_",Prior0spec)))
            return fPDM0help()
        end
        function fIsDrawValid_M0(θ)
            fIDVM0help         =  eval(Meta.parse(string("fIsDrawValid_",Prior0spec)))
            return fIDVM0help(θ)
        end
        function fPriorLogEval_M0(θ)
            fPLEM0help         =  eval(Meta.parse(string("fPriorLogEval_",Prior0spec)))
            return fPLEM0help(θ)
        end


        # Define functions for SMC-MT:

        nP                  = length(vParLabs)

        fLL(θ)            = fLL_M1(θ)
        @everywhere @eval ϕ_Nϕ_M0     = $ϕ_Nϕ_M0
        fLLtilde(θ)       = ϕ_Nϕ_M0 * fLL_M0(θ) + (1-ϕ_Nϕ_M0) * 0

        fPriorDraw()      = fPriorDraw_M1()

        fIsDrawValid(θ)   = minimum(fIsDrawValid_M1(θ))
        fPriorLogEval(θ)  = sum(fPriorLogEval_M1(θ))

    end


    # Setup inputs to SMC algorithm:

    nP              = length(vParLabs)
    showMessages    = 1
    tSettings       = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages)
    @everywhere @eval tSettings     = $tSettings

    setSeed             = 547
    sOutputFilePrefix   = ""
    sOutputPath         = sEstSpecPath


    sSpecFolderM0       = string("dgp",DGPspec,"_m",Model0spec,"_pr",Prior0spec,"_mc",SMC0spec,"/")
    sEstSpecFolderM0    = string( "LTphiLast" , Int(ϕ_Nϕ_M0 * 100) , "/" )

    sEstSpecPathM0      = sMyPath * "Output/" * sSpecFolderM0 * sHelperFolder * sEstSpecFolderM0

    sFileToRead         = "StageLast_Particles.csv"
    sPathToRead         = sEstSpecPathM0 * sFileToRead

    mParticlesFromM0    = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)

    N_M0                = size(mParticlesFromM0,1)
    if N_M0 == N
        mInitParticles      = mParticlesFromM0
    else
        vIndToTake      = rand(Categorical(ones(N)./N),N_M0)
        mInitParticles  = mParticlesFromM0[vIndToTake,:]
    end

    vPIndVariesM1       = fPriorDraw_M1() .!= fPriorDraw_M1()
    vPIndVariesM0       = fPriorDraw_M0() .!= fPriorDraw_M0()
    vPIndvariesOnlyM1   = [vPIndVariesM1[pp] == 1 && vPIndVariesM0[pp] == 0 for pp = 1:nP]
    vPIndDraw           = vPIndvariesOnlyM1


    # Continue interrupted run?

    doesPreviousRunExist = fDoesFileExist(sOutputPath, sOutputFilePrefix* string("StageAll_Stats.csv"))

    if doesPreviousRunExist && continueInterruptedRun == 1

        mContinue_StageAll_Stats    = fGetFile(sOutputPath, sOutputFilePrefix* string("StageAll_Stats.csv"))

        previousRunIncomplete = mContinue_StageAll_Stats[end,1] < 1.0

        if previousRunIncomplete

            nPhiInterrupted             = size(mContinue_StageAll_Stats,1)

            mContinue_WLLP              = fGetFile(sOutputPath * "IntermediateStages/", sOutputFilePrefix* string("Stage",nPhiInterrupted,"_WLLP.csv"))

            mContinue_Particles         = fGetFile(sOutputPath * "IntermediateStages/", sOutputFilePrefix* string("Stage",nPhiInterrupted,"_Particles.csv"))

            # vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,2,setSeed,mInitParticles,vPIndDraw,mContinue_StageAll_Stats,mContinue_WLLP,mContinue_Particles)

            vaMT       = [mInitParticles, vPIndDraw]
            vaContinue = [mContinue_StageAll_Stats, mContinue_WLLP, mContinue_Particles]

            vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,2,setSeed,vaMT,vaContinue)

        end

    end

    if continueInterruptedRun == 0 || doesPreviousRunExist == false

        # Run SMC:

        # vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,2,setSeed,mInitParticles,vPIndDraw)

        vaMT       = [mInitParticles, vPIndDraw]
        vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,2,setSeed,vaMT)

    end


end
