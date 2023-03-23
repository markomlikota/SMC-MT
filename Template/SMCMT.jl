# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# March 2023, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# TEMPLATE SCRIPT TO ESTIMATE SOME MODEL M1 STARTING FROM M0 USING SMC-MT 

# To make it work, adjust lines where comments starting with ### appear



# -------------------------------------------------------------------------------
# OPTIONS:


sTemperingType =            ### "LT" or "MT": 
                            #   LT : likelihood-tempering estimation of M0
                            #   MT : model-tempering estimation of M1 (based on preliminary LT-estimation of M0)

ϕ_Nϕ_M0         =           ### scalar ∈ (0,1], indicating whether LT-estimation of M0 is stopped early (ψ_* in paper;
                            #   fraction of information from M0 which is combined with prior to form proposal)


# Notes:

# - θ is a vector of all parameters appearing in either M0 or M1 (or both).
# - See fSetupSMC.jl for explanations of output of fSMC-function and some more hints on how to specify likelihood and prior functions.
# - fIsDrawValid and fPriorLogEval are specified to return vectors (rather than scalars) because the set of parameters which vary under M1 and M0 can be different
#   (current code allows for additional parameters in M1; extension to have M0-specific parameters, which don't appear in M1, is possible and outlined in paper).
#   In this case, the proposal is composed of the the (tempered) M0-posterior for the common parameters and the prior for the M1-specific parameters, 
#   which is why the above functions need to be able to evaluate and assess the validity of these parameters separately from the rest.
# - To fix some parameters, can either i) drop them from θ, but instead define them to be fixed inside the likelihood function, 
#   or ii) let fPriorDraw_MX return a fixed value in its position, fIsDrawValid give a 1 in its position, and fPriorLogEval give a 0 (log-probability) in its position.
# - For advanced settings, like continuing interrupted runs, or using different SMC-settings (e.g. N) for the two models, see the codes for our applications.
# - For more complicated software-setup, like estimating multiple models (multiple parameter combinations), different priors, different SMC-settings, data samples, etc., 
#   can specify these different settings as "SpecFiles", as in the applications in the paper.
# - The way the functions are written (and the SMC-MT algorithm is discussed in the paper), priors for parameters common to both M0 and M1 need to be the same in the MT-estimation of M1 as in the preliminary LT-estimation of M0.
#   However, presumably, one might be able to accommodate different priors by adjusting the fLL-function on line 289 below as follows:
#       vPIndvariesBoth     = [vPIndVariesM1[pp] == 1 && vPIndVariesM0[pp] == 1 for pp = 1:nP]
#       fLL(θ)              = fLL_M1(θ) + sum( fPriorLogEval_M1(θ)[vPIndvariesBoth] - fPriorLogEval_M0(θ)[vPIndvariesBoth] )



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


nProcs = 8                  # number of workers (parallel computing); could also use just 1

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


# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


@everywhere sMyPath     = pwd() * "/"                           # overall path where software is


include(sMyPath * "/Functions/fSetupSMC.jl")                    # SMC algorithm itself: "fSMC" function
@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")    # for resampling step of SMC
@everywhere include(sMyPath * "/Functions/fHelpersSMC.jl")      # various helping functions for SMC
include(sMyPath * "/Functions/fFolderFileManagement.jl")        # functions to write, read and append csv files


sOutputPathBothModels             =                 ### set path where output shall be stored, ending with "/" ; could use "Output" folder in overall path.
                                                    #   In this path, code will create separate folders for LT-estimation of M0 and MT-estimation of M1: "LTM0", "MTM1".
                                                    #   Could also do that manually.


if sTemperingType == "LT"

    sFolderToCreate = "LTM0"

end

if sTemperingType == "MT"

    sFolderToCreate = "MTM1"

end

try
    mkdir(sOutputPathBothModels * sFolderToCreate)
catch
end
                                                    


# -------------------------------------------------------------------------------

# DEFINE LIKELIHOOD AND PRIOR FUNCTIONS
# -------------------------------------------------------------------------------


mData                       =           ### load data (e.g. (T x N) matrix); could use folder "Data"
@everywhere @eval mData     = $mData    # distribute data to all workers, to be used in likelihood function(s) (implicitly) (and prior; e.g. Minnesota prior for VARs)

@everywhere vParLabs        =           ### define nP x 1 vector of strings (parameter labels), corresponding to scalar-parameters in nP x 1 vector θ



# --- MODEL M1 ---

@everywhere begin                       # each worker needs to know these functions


    function fLL_M1(θ)                      # log-likelihood function for model M1 (model of interest)
        # Input: 
        # - θ :  nP x 1 vector of parameter-values
        # Output: 
        # - LL : log-likelihood (scalar)

        ### ...
        
        return LL
    end

    function fPriorDraw_M1()                # draws θ from prior distribution, model M1
        # Output: 
        # - θ : nP x 1 vector of parameter-values

        ### ...

        return θ
    end

    function fIsDrawValid_M1(θ)             # determines whether a draw θ (from prior or in M-step of SMC) is valid under model M1
        # Input: 
        # - θ :  nP x 1 vector of parameter-values
        # Output: 
        # - vIsValid : nP x 1 vector of 0s and 1s, indicating for each parameter whether it satisfies domain 
        #              (if joint domain is restricted in some particular way, could just return all zeros; 
        #               note that the restrictions placed on parameters common to both M0 and M1 need to be the same under both models)

        ### ...

        return vIsValid
    end

    function fPriorLogEval_M1(θ)            # evaluates prior at θ under model M1
        # Input: 
        # - θ :  nP x 1 vector of parameter-values
        # Output: 
        # - vP : nP x 1 vector of log-probabilities 
        #        (if a set of scalar parameters have a joint distribution (e.g. some matrix Σ distributed as Inverse Wishart), 
        #         put the joint probability in one of the positions and 0s everywhere else)

        ### ...

        return vP
    end

end

# Notes: 
# - Can use matrix "mParaPointers" (output of "fGetParaPointers" in "fHelpersSMC.jl") to map groups of scalar parameters (e.g. Φ,Σ in VAR) to vector θ and vice-versa.
# - Can define additional functions (e.g. in folder "Functions") to be used for likelihood computation (or inside the prior-functions).



# --- MODEL M0 ---

@everywhere begin
        
    function fLL_M0(θ)                      # log-likelihood function for model M0 (approximating model)
        # Input: 
        # - θ :  nP x 1 vector of parameter-values
        # Output: 
        # - LL : log-likelihood (scalar)

        ### ...
        
        return LL
    end

    function fPriorDraw_M0()                # draws θ from prior distribution, model M0
        # Output: 
        # - θ : nP x 1 vector of parameter-values

        ### ...

        return θ
    end

    function fIsDrawValid_M0(θ)             # determines whether a draw θ (from prior or in M-step of SMC) is valid under model M0
        # Input: 
        # - θ :  nP x 1 vector of parameter-values
        # Output: 
        # - vIsValid : nP x 1 vector of 0s and 1s, indicating for each parameter whether it satisfies domain 
        #              (if joint domain is restricted in some particular way, could just return all zeros; 
        #               note that the restrictions placed on parameters common to both M0 and M1 need to be the same under both models)

        ### ...

        return vIsValid
    end

    function fPriorLogEval_M0(θ)            # evaluates prior at θ under model M0
        # Input: 
        # - θ :  nP x 1 vector of parameter-values
        # Output: 
        # - vP : nP x 1 vector of log-probabilities 
        #        (if a set of scalar parameters have a joint distribution (e.g. some matrix Σ distributed as Inverse Wishart), 
        #         put the joint probability in one of the positions and 0s everywhere else)

        ### ...

        return vP
    end

end



# -------------------------------------------------------------------------------

# DEFINE SMC TUNING PARAMETERS
# -------------------------------------------------------------------------------


@everywhere begin
        
    N               = 1000              # number of SMC-particles
    useAdaptive     = 1                 # use adaptive tempering schedule (=1) or fixed TS (=0)
    λ               = 2                 # sets tempering under FTS (ignored under ATS)
    Nϕ              = 200               # sets tempering under FTS (ignored under ATS)
    α               = 0.95              # sets tempering under ATS (ignored under FTS)
    Nmh             = 2                 # number of MH-steps in mutation-step of SMC
    nB              = 1                 # number of blocks for RWMH in mutation steps; this code only supports 1 (no blocking)
    nP              = length(vParLabs)  # number of parameters in θ
    c0              = 0.5               # initial scaling of proposal covariance in mutation-step
    accStar         = 0.25              # target acceptance rate, to iteratively determine scaling of proposal covariance in mutation-step
    Nbar            = N/2               # determines when resampling is triggered
    showMessages    = 1                 # show messages at each SMC-iteration (=1) or not (=0)

    # Notes:
    # - See outline of algorithm in paper and "fSMC" description for more details on what these settings indicate.
    # - Under fixed tempering schedule for LT-estimation of M0, the algorithm terminates exactly when tempering parameter ϕ is equal to ϕ_Nϕ_M0, 
    #   under adaptive tempering schedule, it terminates when ϕ becomes larger than (or equal to) ϕ_Nϕ_M0.


end



# -------------------------------------------------------------------------------

# ESTIMATION
# -------------------------------------------------------------------------------



# LT ESTIMATION OF M0 :


if sTemperingType == "LT"


    # Define functions for SMC-MT:

    @everywhere begin

        fLL(θ)              = fLL_M0(θ)
        fLLtilde(θ)         = 0

        fPriorDraw()        = fPriorDraw_M0()

        fIsDrawValid(θ)     = minimum(fIsDrawValid_M0(θ))
        fPriorLogEval(θ)    = sum(fPriorLogEval_M0(θ))

    end


    # Setup inputs to SMC algorithm:

    tSettings       = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ_M0,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages) # combine SMC-settings into one object; use ϕ_Nϕ = ϕ_Nϕ_M0 here in LT-estimation of M0 
    @everywhere @eval tSettings     = $tSettings                                            # distribute to all workers

    setSeed             =               ### set seed for SMC estimation
    sOutputFilePrefix   = ""            # optional: can let files start with some string to differentiate models e.g.
    sOutputPath         = sOutputPathBothModels * "LTM0/"

    writeFiles          = 1             # see fSMC - description


    # Run SMC:

    vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,writeFiles,setSeed)


end




# MT ESTIMATION OF M1 :


if sTemperingType == "MT"


    # Define functions for SMC-MT:

    @everywhere begin

        fLL(θ)            = fLL_M1(θ)
        @everywhere @eval ϕ_Nϕ_M0     = $ϕ_Nϕ_M0
        fLLtilde(θ)       = ϕ_Nϕ_M0 * fLL_M0(θ) + (1-ϕ_Nϕ_M0) * 0

        fPriorDraw()      = fPriorDraw_M1()

        fIsDrawValid(θ)   = minimum(fIsDrawValid_M1(θ))
        fPriorLogEval(θ)  = sum(fPriorLogEval_M1(θ))

    end


    # Setup inputs to SMC algorithm:


    tSettings       = (N,useAdaptive,λ,Nϕ,1,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages) # combine SMC-settings into one object; use ϕ_Nϕ = 1 here in MT-estimation of M1 
    @everywhere @eval tSettings     = $tSettings                                      # distribute to all workers

    setSeed             =               ### set seed for SMC estimation
    sOutputFilePrefix   = ""            # optional: can let files start with some string to differentiate models e.g.
    sOutputPath         = sOutputPathBothModels * "MTM1/"

    writeFiles          = 1             # see fSMC - description


    # Load posterior particles from LT-estimation of M0:

    sEstSpecPathM0      = sOutputPathBothModels * "LTM0/"

    sFileToRead         = "StageLast_Particles.csv"
    sPathToRead         = sEstSpecPathM0 * sFileToRead

    mParticlesFromM0    = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)

    mInitParticles      = mParticlesFromM0

    vPIndVariesM1       = fPriorDraw_M1() .!= fPriorDraw_M1()
    vPIndVariesM0       = fPriorDraw_M0() .!= fPriorDraw_M0()
    vPIndvariesOnlyM1   = [vPIndVariesM1[pp] == 1 && vPIndVariesM0[pp] == 0 for pp = 1:nP]
    vPIndDraw           = vPIndvariesOnlyM1


    # Run SMC:
    
    vaMT       = [mInitParticles, vPIndDraw]
    vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,writeFiles,setSeed,vaMT)


end
