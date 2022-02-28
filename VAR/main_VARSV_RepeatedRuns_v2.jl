# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Main script to analyze estimation of VAR-SV using SMC-LT/MT over many runs.



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec         = "1"     # data from which DGP to take?
SMCspec         = "1"     # which values for SMC tuning parameters to take?

Model1spec      = "SVx8"      # which model (LL) to estimate?
Prior1spec      = "SVx81"     # which prior to take?

Model0spec      = "HM"            # which model is M0?
Prior0spec      = "HM1"           # which prior to take for M0?

vPhiLast        = [0.2, 0.4, 0.6, 0.8, 1.0]

# Remaining options (see main_VARSV_OneRun.jl) specified below inside function fRunSMC

nRun            = 1

sRunFolderName  = "xrun"

seed            = 2 # 0: no seed, otherwise any positive integer



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
using Roots
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



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


# @everywhere cd()
# @everywhere sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareVAR/")
@everywhere sMyPath            = pwd() * "/"


@everywhere include(sMyPath * "/Functions/fHelpersSMC.jl")
include(sMyPath * "/Functions/fVAR.jl")
include(sMyPath * "/Functions/fSetupSMC.jl")
@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")
include(sMyPath * "/Functions/fFolderFileManagement.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# -- Load Data, Model, Prior, SMC-info -----------------------------------------


# Load data:

sFileToRead         = "Data_" * string("dgp",DGPspec) * ".csv"
sPathToRead         = sMyPath * "Data/" * sFileToRead

# # mData               = CSV.read(sPathToRead;header=1) # DOESN'T WORK YET ...
# mData               = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)
mData               = fMyCSVREAD(sPathToRead)
@everywhere @eval mData     = $mData


# Load Model (Likelihood) and prior:

@everywhere @eval Model1spec     = $Model1spec
@everywhere @eval Prior1spec     = $Prior1spec

@everywhere include(sMyPath * "SpecFiles/" * string("script_Modelspec",Model1spec,".jl") )
@everywhere include(sMyPath * "SpecFiles/" * string("script_Priorspec",Prior1spec,".jl") )

@everywhere @eval Model0spec     = $Model0spec
@everywhere @eval Prior0spec     = $Prior0spec

@everywhere include(sMyPath * "SpecFiles/" * string("script_Modelspec",Model0spec,".jl") )
@everywhere include(sMyPath * "SpecFiles/" * string("script_Priorspec",Prior0spec,".jl") )


# Load SMC tuning parameters:

include(sMyPath * "SpecFiles/" * string("script_SMCspec",SMCspec,".jl") )



# -- Define Function to perform SMC estimation repeatedly ----------------------

function fRunSMC(DGPspec,M1spec,PM1spec,SMCspec,sTemperingType,ϕ_Nϕ,createNewRunFolder,sRunFolderName,setSeed,M0spec="HM",PM0spec="HM1",ϕ_Nϕ_M0=1.0)


    @everywhere @eval M1spec     = $M1spec
    @everywhere @eval PM1spec     = $PM1spec
    @everywhere @eval M0spec     = $M0spec
    @everywhere @eval PM0spec     = $PM0spec


    # -- Set up directory for estimation output ------------------------------------


    sSpecFolder             = string("dgp",DGPspec,"_m",M1spec,"_pr",PM1spec,"_mc",SMCspec,"/")
    sSpecPath               = sMyPath * "Output/" * sSpecFolder
    try
        mkdir(sSpecPath)
    catch
    end

    if createNewRunFolder == 1
        sHelperFolder           = fAppendFolder(sSpecPath,sRunFolderName) * "/"
    else
        sHelperFolder           = fFindNewestFolder(sSpecPath,sRunFolderName)[1] * "/"
    end
    sHelperPath             = sSpecPath * sHelperFolder
    try
        mkdir(sHelperPath)
    catch
    end

    if sTemperingType == "LT"
        sEstSpecFolder          = string( "LTphiLast" , Int(ϕ_Nϕ * 100) , "/" )
    elseif sTemperingType == "MT"
        sEstSpecFolder          = string( "MTphiLast" , Int(ϕ_Nϕ_M0 * 100) , "_m", M0spec, "_pr", PM0spec, "/" )
    end
    sEstSpecPath            = sHelperPath * sEstSpecFolder
    try
        mkdir(sEstSpecPath)
    catch
    end



    # -- Estimation ----------------------------------------------------------------


    if sTemperingType == "LT"

        # Define LL, LL-tilde and prior functions :

        @everywhere begin

            #Need to use anonymous functions because they are changed every time Model or prior spec files are changed.

            vParLabs        = eval(Meta.parse(string("vParLabs_",M1spec)))
            function fLL(θ)
                fLLhelp     = eval(Meta.parse(string("fLL_",M1spec)))
                return fLLhelp(θ)
            end
            fLLtilde(θ)     = 0
            function fPriorDraw()
                fPDhelp     = eval(Meta.parse(string("fPriorDraw_",M1spec)))
                return fPDhelp()
            end
            function fIsDrawValid(θ)
                fIDVhelp         =  eval(Meta.parse(string("fIsDrawValid_",M1spec)))
                return minimum(fIDVhelp(θ))
            end
            function fPriorLogEval(θ)
                fPLEhelp         =  eval(Meta.parse(string("fPriorLogEval_",M1spec)))
                return sum(fPLEhelp(θ))
            end

        end


        # Setup inputs to SMC algorithm:

        nP              = length(vParLabs)
        tSettings       = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages)
        @everywhere @eval tSettings     = $tSettings

        sOutputFilePrefix   = ""
        sOutputPath         = sEstSpecPath


        # Run SMC:

        fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed)


    elseif sTemperingType == "MT"

        # Define LL, LL-tilde and prior functions :
        # (written in most general way, but stuff simplifies in our VAR-SV application because:
        # a) θ^0 = ∅ (no parameters unique to M0), and b) we take same prior for θ^v under both models)

        @everywhere begin

            # First "translate" model functions to M1 and M0:

            vParLabs_M1        = eval(Meta.parse(string("vParLabs_",M1spec)))

            function fLL_M1(ϑM1)
                fLLM1help     = eval(Meta.parse(string("fLL_",M1spec)))
                return fLLM1help(ϑM1)
            end
            function fPriorDraw_M1()
                fPDM1help     = eval(Meta.parse(string("fPriorDraw_",M1spec)))
                return fPDM1help()
            end
            function fIsDrawValid_M1(ϑM1)
                fIDVM1help         =  eval(Meta.parse(string("fIsDrawValid_",M1spec)))
                return fIDVM1help(ϑM1)
            end
            function fPriorLogEval_M1(ϑM1)
                fPLEM1help         =  eval(Meta.parse(string("fPriorLogEval_",M1spec)))
                return fPLEM1help(ϑM1)
            end


            vParLabs_M0        = eval(Meta.parse(string("vParLabs_",M0spec)))

            function fLL_M0(ϑM0)
                fLLM0help     = eval(Meta.parse(string("fLL_",M0spec)))
                return fLLM0help(ϑM0)
            end
            function fPriorDraw_M0()
                fPDM0help     = eval(Meta.parse(string("fPriorDraw_",M0spec)))
                return fPDM0help()
            end
            function fIsDrawValid_M0(ϑM0)
                fIDVM0help         =  eval(Meta.parse(string("fIsDrawValid_",M0spec)))
                return fIDVM0help(ϑM0)
            end
            function fPriorLogEval_M0(ϑM0)
                fPLEM0help         =  eval(Meta.parse(string("fPriorLogEval_",M0spec)))
                return fPLEM0help(ϑM0)
            end


            # Define functions for SMC-MT:

            vParLabs        = vParLabs_M1

            nP              = length(vParLabs)
            nPM1            = length(vParLabs_M1)
            nPM0            = length(vParLabs_M0)

            vPInd_VariesM1  = [trues(nP)...]
            vPInd_VariesM0  = [trues(nPM0)...,falses(nPM1-nPM0)...]

            vPInd_VariesOnlyM0 = [vPInd_VariesM1[pp] .== false && vPInd_VariesM0[pp] .== true for pp = 1:nP]        # go from θ to θ^0
            vPInd_VariesOnlyM1 = [vPInd_VariesM1[pp] .== true && vPInd_VariesM0[pp] .== false for pp = 1:nP]        # go from θ to θ^1
            vPInd_VariesBoth   = [vPInd_VariesM1[pp] .== true && vPInd_VariesM0[pp] .== true for pp = 1:nP]         # go from θ to θ^v

            vPIndM0_VariesOnlyM0 = vPInd_VariesOnlyM0[vPInd_VariesM0]     # go from ϑM0 to θ^0
            vPIndM0_VariesBoth = vPInd_VariesBoth[vPInd_VariesM0]         # go from ϑM0 to θ^v
            vPIndM1_VariesBoth = vPInd_VariesBoth[vPInd_VariesM1]         # go from ϑM1 to θ^v


            fLL(θ)            = fLL_M1(θ[vPInd_VariesM1]) #+
                                # sum( fPriorLogEval_M1(θ[vPInd_VariesM1])[vPIndM1_VariesBoth] -
                                #      fPriorLogEval_M0(θ[vPInd_VariesM0])[vPIndM0_VariesBoth] )
            @eval ϕ_Nϕ_M0     = $ϕ_Nϕ_M0
            fLLtilde(θ)       = ϕ_Nϕ_M0 * fLL_M0(θ[vPInd_VariesM0]) + (1-ϕ_Nϕ_M0) * 0

            fPriorDraw()      = [fPriorDraw_M1()..., fPriorDraw_M0()[vPIndM0_VariesOnlyM0]...]

            fIsDrawValid(θ)    = minimum( [fIsDrawValid_M1(θ[vPInd_VariesM1])...,
                                          fIsDrawValid_M0(θ[vPInd_VariesM0])[vPIndM0_VariesOnlyM0]...] )
            fPriorLogEval(θ)  = sum( [fPriorLogEval_M1(θ[vPInd_VariesM1])...,
                                      fPriorLogEval_M0(θ[vPInd_VariesM0])[vPIndM0_VariesOnlyM0]...] )

        end


        # Setup inputs to SMC algorithm:

        nP              = length(vParLabs)
        tSettings       = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages)
        @everywhere @eval tSettings     = $tSettings

        sOutputFilePrefix   = ""
        sOutputPath         = sEstSpecPath


        sSpecFolderM0       = string("dgp",DGPspec,"_m",M0spec,"_pr",PM0spec,"_mc",SMCspec,"/")
        sEstSpecFolderM0    = string( "LTphiLast" , Int(ϕ_Nϕ_M0 * 100) , "/" )

        sEstSpecPathM0      = sMyPath * "Output/" * sSpecFolderM0 * sHelperFolder * sEstSpecFolderM0

        sFileToRead         = "StageLast_Particles.csv"
        sPathToRead         = sEstSpecPathM0 * sFileToRead

        # mParticlesFromM0    = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)
        mParticlesFromM0    = fMyCSVREAD(sPathToRead)

        mInitParticles      = zeros(N,nP)
        mInitParticles[:,vPInd_VariesM0] = mParticlesFromM0

        vPIndDraw          = vPInd_VariesOnlyM1


        # Run SMC:

        fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed,mInitParticles,vPIndDraw)


    end


end



# -- Execute it repeatedly and for different estimation specifications ---------


for rr = 1:nRun

    println("*************************")
    println("*************************")
    println(string("     RUN     ",rr,"     "))
    println("*************************")
    println("*************************")


    # Estimate M1 (=SV) using LT:

    println("")
    println(string("Run ", rr, ": M1 - LT "))
    println("")

    Model1specHere      = deepcopy(Model1spec)
    Prior1specHere      = deepcopy(Prior1spec)
    sTemperingType      = "LT"
    ϕ_Nϕ                = 1.0

    createNewRunFolder = 1
    setSeed             =  seed*rr*1

    @time fRunSMC(DGPspec,Model1specHere,Prior1specHere,SMCspec,sTemperingType,ϕ_Nϕ,createNewRunFolder,sRunFolderName,setSeed)


    for pphi = 1:length(vPhiLast)

        ϕ_Nϕ_M0         = vPhiLast[pphi]


        # Estimate M0 (=HM) using LT:

        println("")
        println(string("Run ", rr, ": M0 - LT - ϕ_Nϕ = ",ϕ_Nϕ_M0))
        println("")

        Model1specHere  = deepcopy(Model0spec)
        Prior1specHere  = deepcopy(Prior0spec)
        sTemperingType  = "LT"
        ϕ_Nϕ            = deepcopy(ϕ_Nϕ_M0)

        pphi == 1 ? createNewRunFolder = 1 : createNewRunFolder = 0

        setSeed             = seed*rr*(pphi+1)*3

        @time fRunSMC(DGPspec,Model1specHere,Prior1specHere,SMCspec,sTemperingType,ϕ_Nϕ,createNewRunFolder,sRunFolderName,setSeed)


        # Estimate M1 (=SV) using MT:

        println("")
        println(string("Run ", rr, ": M1 - MT - ϕ_Nϕ_M0 = ",ϕ_Nϕ_M0))
        println("")

        Model1specHere      = deepcopy(Model1spec)
        Prior1specHere      = deepcopy(Prior1spec)
        sTemperingType      = "MT"
        ϕ_Nϕ                = 1.0
        Model0specHere      = deepcopy(Model0spec)
        Prior0specHere      = deepcopy(Prior0spec)

        createNewRunFolder = 0

        setSeed             = seed*rr*(pphi+1)*5

        @time fRunSMC(DGPspec,Model1specHere,Prior1specHere,SMCspec,sTemperingType,ϕ_Nϕ,createNewRunFolder,sRunFolderName,setSeed,Model0specHere,Prior0specHere,ϕ_Nϕ_M0)

    end

end
