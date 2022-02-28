# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Given target distribution (& SMC settings), runs SMC for many different proposal distributions, each many (nRun) times.

# Code for bivariate case only.

# Note: only target-statistics (mean, std, etc.) for first parameter are computed after each SMC run.

# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


Targetspec          = "B0"
SMCspec             = "1"

nRun                = 1 # for each proposal density considered

vMeans01             = -3.0:1.0:0.0 #-3:0.5:0
vStds01              = 1.0:1.0:2.0 #0.2:0.4:2

vMeans02            = [-1,0]
vStds02             = [1,2]
vRhos0              = [-0.5,0] #[-0.5,0,0.5]

setOfProposals      = 2 # either 0 or integer,
# 0: create and store output in new folder "SetOfProposalsX+1" (appends list of folders),
# 2 (e.g.): store results in output folder "SetOfProposals2" (e.g. when adding runs to existing folder)(Note: vMeans01 and vStds01 are ignored then; taken from existing folder instead)
# KEEP AT 0 BY DEFAULT!!! (not to add add results to wrong folder)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


# nProcs = 8
#
# using Distributed
# addprocs(nProcs-length(procs()))
#
#
# using Plots
# pyplot()
# using Roots
# using StatsBase
# using KernelDensity
# using DataFrames
# using CSV
# using StatsFuns
# using DelimitedFiles
#
# @everywhere using ParallelDataTransfer
# @everywhere using PositiveFactorizations
# @everywhere using Random
# @everywhere using SharedArrays
# @everywhere using Distributions
# @everywhere using LinearAlgebra



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


@everywhere cd()
@everywhere sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareSimpleExample/")


include(sMyPath * "/Functions/fSetupSMC.jl")
@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")
@everywhere include(sMyPath * "/Functions/fHelpers.jl")
include(sMyPath * "/Functions/fFolderFileManagement.jl")
include(sMyPath * "/Functions/fMyPlot.jl")
@everywhere include(sMyPath * "/Functions/fGetProposalObjects.jl")
include(sMyPath * "/Functions/fGetDistanceMeasures.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------



# -- Set up directory for estimation output ------------------------------------


sSpecFolder             = string("t",Targetspec,"_smc",SMCspec,"/")
sSpecPath               = sMyPath * "Output/" * sSpecFolder

if setOfProposals == 0
    sHelperFolder        = fAppendFolder(sSpecPath,"SetOfProposals") * "/"
else
    sHelperFolder        = string("SetOfProposals",setOfProposals) * "/"
end
sHelperPath     = sSpecPath * sHelperFolder

sEstSpecFolder          = ""
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



# -- Estimation -----------------------------------------------------------------


# Set up target density and SMC settings:

@everywhere @eval Targetspec  = $Targetspec

@everywhere include(sMyPath * "/SpecFiles/" * string("script_Targetspec",Targetspec,".jl") )

include(sMyPath * "/SpecFiles/" * string("script_SMCspec",SMCspec,".jl") )




# Setup output folder and filename prefix:

sOutputPath       = sEstSpecPath
sOutputFilePrefix = string("t",Targetspec,"_smc",SMCspec,"_")


# Write/Read file with proposal dist. parameters (μ0,σ0) and corresp. distance measures :

cd(sOutputPath)

if setOfProposals == 0

    sSaveName           = sOutputFilePrefix * "parameters.csv"
    vColLabs1           = string.("p",1:5)
    vColLabs2           = ["mean1","sig1","mean2","sig2","rho"]
    vColLabs            = [vColLabs1...,vColLabs2...,"DM_areaprod","DM_KL","DM_areamin","DM_IS"]
    CSV.write(sOutputPath * sSaveName,  DataFrame(permutedims(vColLabs)), header=false)

    nMeans01            = length(vMeans01)
    nStds01             = length(vStds01)
    nMeans02            = length(vMeans02)
    nStds02             = length(vStds02)
    nRhos0              = length(vRhos0)
    mSpaceInd           = fGetCartProd([1:nMeans01,1:nStds01,1:nMeans02,1:nStds02,1:nRhos0])[1]
    mSpaceParams        = fGetCartProd([vMeans01,vStds01,vMeans02,vStds02,vRhos0])[1]


    mParamInfo          = [mSpaceInd mSpaceParams]

    nProps              = size(mParamInfo,1)

    mDistanceMeasures   = zeros(nProps,4)
    for kk = 1:nProps
        μ1g,σ1g,μ2g,σ2g,ρg = mSpaceParams[kk,:]
        covg     = ρg * σ1g * σ2g
        mΣg      = [σ1g^2   covg    ;
                    covg    σ2g^2   ]
        vμg      = [μ1g,μ2g]
        mDistanceMeasures[kk,:] = fGetDistanceMeasures(vμf,mΣf,vμg,mΣg)
    end

    CSV.write(sOutputPath * sSaveName,  DataFrame([mParamInfo mDistanceMeasures]), append=true)

else

    sReadName            = sOutputFilePrefix * "parameters.csv"
    # dfParams             = CSV.read(sOutputPath * sReadName, DataFrame;header=1)
    # mParams              = Matrix{Float64}(dfParams)
    mParams              = readdlm(sOutputPath * sReadName, ',', Float64, '\n'; skipstart=1)

    mParamInfo           = mParams[:,1:10]
    mDistanceMeasures    = mParams[:,11:end]
    nProps              = size(mParamInfo,1)

    vMeans01              = unique(mParamInfo[:,6])
    vStds01               = unique(mParamInfo[:,7])
    vMeans02             = unique(mParamInfo[:,8])
    vStds02              = unique(mParamInfo[:,9])
    vRhos0               = unique(mParamInfo[:,10])

    nMeans01             = length(vMeans01)
    nStds01              = length(vStds01)
    nMeans02            = length(vMeans02)
    nStds02             = length(vStds02)
    nRhos0              = length(vRhos0)

end



# Run estimation:

vStatsNames         = ["mean","std","skewness","kurtosis","p2.5","p5","p95","p97.5","'MDD'","nStages","time"]

for rr = 1:nRun
for kk = 1:nProps

    mm1,ss1,mm2,ss2,pp = Int.(mParamInfo[kk,1:5])
    μ1g,σ1g,μ2g,σ2g,ρg = mParamInfo[kk,6:10]

    # mm <= 5 ? continue : nothing
    # mm == 6 && ss < 2 ? continue : nothing
    # mm == 6 && ss==2 && rr < 68 ? continue : nothing


    println("*************************")
    println(string(" Mean1  ",mm1,"/",nMeans01,", Std1  ",ss1,"/",nStds01,", Mean2  ",mm2,"/",nMeans02,", Std2  ",ss2,"/",nStds02,", Rho  ",pp,"/",nRhos0,", Run  ",rr,"/",nRun))
    println("*************************")


    # Setup filename for proposal and write mean and std to one file and tempering schedule to another:

    sProposalFileName   = sOutputFilePrefix * "p_" * fCombinevStrings(Int.(mParamInfo[kk,1:5]))
    vFileNames          = readdir(sOutputPath)

    sSaveName           = sProposalFileName * "_statstics.csv"
    vInd                = findall(occursin.(sSaveName,vFileNames) .== 1)
    doesFileExist       = length(vInd) > 0
    if doesFileExist == false # only if file does not exist, write it anew; this allows to add more runs easily to given model (proposal) specifications
      CSV.write(sOutputPath * sSaveName,  DataFrame(permutedims(vStatsNames)), header=false)
    end

    sSaveName           = sProposalFileName * "_temperingschedule.csv"
    vInd                = findall(occursin.(sSaveName,vFileNames) .== 1)
    doesFileExist       = length(vInd) > 0
    if doesFileExist == false # only if file does not exist, write it anew; this allows to add more runs easily to given model (proposal) specifications
      CSV.write(sOutputPath * sSaveName,  DataFrame(permutedims(["stage1","stage2"])), header=false)
    end


    # Define functions/objects for SMC:

    @eval @everywhere μ1g=$μ1g
    @eval @everywhere σ1g=$σ1g
    @eval @everywhere μ2g=$μ2g
    @eval @everywhere σ2g=$σ2g
    @eval @everywhere ρg=$ρg

    @everywhere begin

      fDrawProposal, fEvalLogProposal, dProposal = fGetProposalObjects_d2(μ1g,σ1g,μ2g,σ2g,ρg)

      fLLtilde(θ)         = fEvalLogProposal(θ)
      fLL(θ)              = fEvalLogTarget(θ)
      #fPriorDraw()        = fPriorDrawM1() #not needed here; (all) initial particles supplied
      fIsDrawValid(θ)     = 1
      fPriorLogEval(θ)    = 0

    end

    vParLabs  = string.("p",1:nP)
    @eval @everywhere vParLabs = $vParLabs


    # Draw initial particles from proposal distribution:

    mInitParticles        = zeros(N,nP)
    [mInitParticles[nn,:] = fDrawProposal() for nn = 1:N]

    vPIndDraw               = falses(nP)


    # Run SMC to obtain target distribution:

    setSeed               = 0

    vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed,mInitParticles,vPIndDraw,0)

    if kk==1 && rr==1 #run twice for Julia to compile (affects runtime only; first time is slower)
        vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed,mInitParticles,vPIndDraw,0)
    end


    # Compute and write statistics of interest for this run:

    vStatsHere      = zeros(length(vStatsNames),1)

    vPosteriorFirstParam = vmParticles[end][:,1]

    vStatsHere[1]   = mean(vPosteriorFirstParam)
    vStatsHere[2]   = std(vPosteriorFirstParam)
    vStatsHere[3]   = skewness(vPosteriorFirstParam)
    vStatsHere[4]   = kurtosis(vPosteriorFirstParam)


    vParticlesEndSorted = sort(vPosteriorFirstParam)
    vStatsHere[5]   = vParticlesEndSorted[Int(N*0.025)]
    vStatsHere[6]   = vParticlesEndSorted[Int(N*0.05)]
    vStatsHere[7]   = vParticlesEndSorted[Int(N*0.95)]
    vStatsHere[8]   = vParticlesEndSorted[Int(N*0.975)]

    vStatsHere[9]   = sum(vLogMDD)
    vStatsHere[10]  = length(vmParticles)-1
    vStatsHere[11]  = timeTotal

    sSaveName       = sProposalFileName * "_statstics.csv"
    CSV.write(sOutputPath * sSaveName,  DataFrame(vStatsHere'), append=true)

    sSaveName       = sProposalFileName * "_temperingschedule.csv"
    CSV.write(sOutputPath * sSaveName,  DataFrame(vϕ'), append=true)


end
end
