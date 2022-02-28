# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Feb 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Reads in posterior draws from one specification and computes LL at these parameter values under different models



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


# Where to take posterior draw from:

DGPspec         = "3"
Model1spec      = "L2"
Prior1spec      = "L22"
SMC1spec        = "03"

# sTemperingType  will be "LT" & correspondingly ϕ_Nϕ = 1.0

Model0spec      = "L1"



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
using KernelDensity
using DataFrames
using CSV
using DelimitedFiles
using Statistics

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


@everywhere cd()
@everywhere sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/RBC/")
# sMyPath            = pwd() * "/"


include(sMyPath * "Functions/fHelpersSMC.jl")
include(sMyPath * "Functions/fFolderFileManagement.jl")
include(sMyPath * "Functions/fFunctionsOutputAnalysis.jl")
@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


sTemperingType  = "LT"
ϕ_Nϕ            = 1.0


# -- Preliminary ---------------------------------------------------------------


sSpecPath       = sMyPath * "Output/" * fGetSpecFolder(DGPspec,Model1spec,Prior1spec,SMC1spec)

sFilePath       = sSpecPath * fGetEstSpecFolderLT(ϕ_Nϕ)



# -- Posterior Evolution & Prior-Posterior Plots -------------------------------


# Object with Posterior Particles:

sFileSuffix             = "StageLast_Particles.csv"

mParticlesStageLast     = fGetFile(sFilePath,sFileSuffix)


# Load data:

sFileToRead             = "Data_" * string("dgp",DGPspec) * ".csv"
sPathToRead             = sMyPath * "Data/" * sFileToRead

# mData                 = CSV.read(sPathToRead;header=1) # DOESN'T WORK YET ...
mData                   = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)


# Load models (LL-functions):

@everywhere @eval Model1spec     = $Model1spec

@everywhere include(sMyPath * "SpecFiles/" * string("script_Modelspec",Model1spec,".jl") )


@everywhere @eval Model0spec     = $Model0spec

@everywhere include(sMyPath * "SpecFiles/" * string("script_Modelspec",Model0spec,".jl") )



# Compute runtimes:

@eval @everywhere mParticlesStageLast    = $mParticlesStageLast
@eval @everywhere mData                  = $mData

N           = size(mParticlesStageLast,1)
mTimes      = SharedArray{Float64}((N,2))

@sync @distributed for nn = 1:N
    println(nn)
    θ               = mParticlesStageLast[nn,:]
    mTimes[nn,1]    = @elapsed fLL_L1(θ)
    mTimes[nn,2]    = @elapsed fLL_L2(θ)
end

vTimes      = mean(mTimes,dims=1)
timeM0      = vTimes[1]
timeM1      = vTimes[2]

println("")
println("M0    = ",timeM0,"s")
println("M1    = ",timeM1,"s")
println("Ratio = ",timeM1/timeM0,"s")
