# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Given target distribution (& SMC settings), runs SMC for many different proposal distributions, each many (nRun) times.

# Code for univariate case only.

# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


Targetspec          = "A0" # A0 = standard normal
SMCspec             = "1"

nRun                = 100 # number of repetitions for each proposal density considered

vMeans0             = -3:0.5:0
vStds0              = 0.2:0.2:2

setOfProposals      = 0 # either 0 or integer,
# 0: create and store output in new folder "SetOfProposalsX+1" (appends list of folders),
# 2 (e.g.): store results in output folder "SetOfProposals2" (e.g. when adding runs to existing folder)(Note: vMeans0 and vStds0 are ignored then; taken from existing folder instead)
# KEEP AT 0 BY DEFAULT!!! (not to add add results to wrong folder)



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
using KernelDensity
using DataFrames
using CSV
using StatsFuns
using DelimitedFiles

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
# @everywhere sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareSimpleExample/")
@everywhere sMyPath            = string(pwd(),"/")


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

try
    mkdir(sSpecPath)
catch
end

if setOfProposals == 0
    sHelperFolder       = fAppendFolder(sSpecPath,"SetOfProposals") * "/"
else
    sHelperFolder       = string("SetOfProposals",setOfProposals) * "/"
end
sHelperPath             = sSpecPath * sHelperFolder

sEstSpecFolder          = ""
sEstSpecPath            = sHelperPath * sEstSpecFolder

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

# cd(sOutputPath)

if setOfProposals == 0

    sSaveName            =  sOutputFilePrefix * "parameters.csv"
    CSV.write(sOutputPath * sSaveName,  DataFrame(permutedims(["pNumber1","pNumber2","mean","std","DM_areaprod","DM_KL","DM_areamin","DM_IS"])), header=false)

    nMeans0             = length(vMeans0)
    nStds0              = length(vStds0)
    mParamInfo          = [repeat(1:nMeans0,nStds0,1) repeat(1:nStds0, inner=(nMeans0,1)) repeat(vMeans0,nStds0,1) repeat(vStds0,inner=(nMeans0,1)) ]

    mDistanceMeasures   = zeros(nMeans0*nStds0,4)
    for kk = 1:nMeans0*nStds0
        mm = Int(mParamInfo[kk,1])
        ss = Int(mParamInfo[kk,2])
        μg = vMeans0[mm]
        σg = vStds0[ss]
        dProposal = Normal(μg,σg)
        mDistanceMeasures[kk,:] = fGetDistanceMeasures(μf,σf^2,μg,σg^2)
    end
    # unique( mDistanceMeasures[mParamInfo[:,4].^2 .< σf^2 / 2,4]
    # sum( isnan.( mDistanceMeasures[mParamInfo[:,4].^2 .> σf^2 / 2,4] ) )

    CSV.write(sOutputPath * sSaveName,  DataFrame([mParamInfo mDistanceMeasures]), append=true)

else

    sReadName            = sOutputFilePrefix * "parameters.csv"
    # dfParams             = CSV.read(sOutputPath * sReadName, DataFrame;header=1)
    # mParams              = Matrix{Float64}(dfParams)
    mParams              = readdlm(sOutputPath * sReadName, ',', Float64, '\n'; skipstart=1)

    mParamInfo           = mParams[:,1:4]
    mDistanceMeasures    = mParams[:,5:end]

    vStds0               = unique(mParamInfo[:,4])
    nStds0               = length(vStds0)
    vMeans0              = unique(mParamInfo[:,3])
    nMeans0              = length(vMeans0)

end



# Run estimation:

vStatsNames         = ["mean","std","skewness","kurtosis","p2.5","p5","p95","p97.5","'MDD'","nStages","time"]

for rr = 1:nRun
for mm = 1:nMeans0, ss = 1:nStds0

    # mm <= 5 ? continue : nothing
    # mm == 6 && ss < 2 ? continue : nothing
    # mm == 6 && ss==2 && rr < 68 ? continue : nothing


    println("*************************")
    println(string(" Mean ",mm,"/",nMeans0,", Std ",ss,"/",nStds0,", Run ",rr,"/",nRun))
    println("*************************")


    # Setup filename for proposal and write mean and std to one file and tempering schedule to another:

    sProposalFileName   = string(sOutputFilePrefix,"p_m",mm,"_s",ss)
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

    μg = vMeans0[mm]
    σg = vStds0[ss]

    @eval @everywhere μg=$μg
    @eval @everywhere σg=$σg

    @everywhere begin

      fDrawProposal, fEvalLogProposal, dProposal = fGetProposalObjects(μg,σg)

      fLLtilde(θ)         = fEvalLogProposal(θ)
      fLL(θ)              = fEvalLogTarget(θ)
      #fPriorDraw()        = fPriorDrawM1() #not needed here; (all) initial particles supplied
      fIsDrawValid(θ)     = 1
      fPriorLogEval(θ)    = 0

    end

    vParLabs  = string.("p",1:nP)
    @eval @everywhere vParLabs = $vParLabs


    # Draw initial particles from proposal distribution:

    if nP == 1
      mInitParticles        = [fDrawProposal() for n = 1:N]
      mInitParticles        = reshape(mInitParticles,(N,nP))
    elseif nP == 2
      mInitParticles        = zeros(N,nP)
      [mInitParticles[nn,:] = fDrawProposal() for nn = 1:N]
    end

    vPIndDraw               = falses(nP)


    # Run SMC to obtain target distribution:

    setSeed               = 0

    vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed,mInitParticles,vPIndDraw,0)

    if mm==1 && ss==1 && rr==1 #run twice for Julia to compile (affects runtime only; first time is slower)
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
