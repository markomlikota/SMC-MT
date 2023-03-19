# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# July 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Script to calibrate number of BSPF particles for LL-eval. under mL2



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec             = "1"       # data from which DGP to take?
useDataModelspec    = "L2"      # use "L2" to calibrate BSPF, use "L1" to check whether in this case LL under L2 is lower than under L1

Mbspf               = 25000

addLLunderL1        = 1 # (to plot)


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


@everywhere cd()
@everywhere sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/DSGE-ABS/")
# @everywhere sMyPath            = pwd() * "/"


@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")
include(sMyPath * "/Functions/fMyPlot.jl")
include(sMyPath * "/Functions/fFolderFileManagement.jl")


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


sFolderForOutput    = sMyPath * "Output/Checks/"


# -- Load Data- and Model-info -----------------------------------------


# Load data:

sFileToRead         = "Data_" * string("dgp",DGPspec) * ".csv"
sPathToRead         = sMyPath * "Data/" * sFileToRead

# mData               = CSV.read(sPathToRead;header=1) # DOESN'T WORK YET ...
mData               = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)
@eval @everywhere mData     = $mData


# Load Model :

@everywhere include(sMyPath * "SpecFiles/" * string("script_ModelspecL2.jl") )


# to get θ0:

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )



# -- Investigate BSPF ----------------------------------------------------------


# Mbspf = 500
# loglik, timeSolution, timePF = fLL_L2(θ0,Mbspf)


mData =  fGetFile(sMyPath * "Data/","Data_" * string("dgp",DGPspec) * string("_m",useDataModelspec) * ".csv")

@eval @everywhere Mbspf = $Mbspf
@eval @everywhere θ0    = $θ0
@eval @everywhere mData = $mData

M      = 100
vDraws = SharedArray{Float64}((M))
mTimes = SharedArray{Float64}((M,2))

@sync @distributed for mm = 1:M
    println(mm)
    vDraws[mm], mTimes[mm,1], mTimes[mm,2] = fLL_L2(θ0,Mbspf)
end


vTimes      = mean(mTimes,dims=1)
timeSol     = vTimes[1]
timeBSPF    = vTimes[2]
ppp         = histogram(vDraws,fillalpha=0.5,label="mL2")
title!(string("LL at θ0: BSPF (M=",Mbspf," particles) \n solution = ",round(timeSol;digits=6),"s, BSPF = ",round(timeBSPF;digits=6),"s"))
tYLimsHere = ylims(ppp)

if addLLunderL1 == 1
    include(sMyPath * "SpecFiles/" * string("script_ModelspecL1.jl") )
    timeKF = @elapsed llKF = fLL_L1(θ0)
    timeKF = @elapsed llKF = fLL_L1(θ0) #twice, to compile
    plot!(ones(2,1)*llKF,[tYLimsHere[1],tYLimsHere[2]],line=(:red,0.9,2,:line),label="mL1")
end

display(ppp)

sFileForOutput         = string("plot_mL2_BSPFcheck_M",Mbspf,"_dgp",DGPspec,"_m",useDataModelspec,".png")
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)



timeM0 = timeKF
timeM1 = timeSol + timeBSPF

timeM0 / timeM1
