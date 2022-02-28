# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Script to check whether BSPF works by comparing Models "L1" (linearized RBC, LL evaluated with Kalman Filter) and "L1bspf"



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec         = "1"     # data from which DGP to take?

Mbspf       = 200



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
@everywhere using Roots



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


@everywhere cd()
@everywhere sMyPath            = string(pwd(),"/Dropbox/ModelTempering/Software/RBC/")
# @everywhere sMyPath            = pwd() * "/"


@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")


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
@everywhere @eval mData     = $mData


# Load Models (Likelihood evaluations) :

@everywhere include(sMyPath * "SpecFiles/" * string("script_ModelspecL1.jl") )
@everywhere include(sMyPath * "SpecFiles/" * string("script_ModelspecL1bspf.jl") )



# -- Compare LLs -----------------------------------------


# to get θ0:

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )

M = 1000
vDraws = SharedArray{Float64}((M))
vTimes = SharedArray{Float64}((M))
@everywhere @eval Mbspf = $Mbspf
@everywhere @eval θ0 = $θ0

@sync @distributed for mm = 1:M
    println(mm)
    vTimes[mm] = @elapsed vDraws[mm] = fLL_L1bspf(θ0,Mbspf)
end

timeKF = @elapsed llKF = fLL_L1(θ0)
timeKF = @elapsed llKF = fLL_L1(θ0) #twice, to compile

timeBSPF = mean(vTimes)
ppp = histogram(vDraws,fillalpha=0.5,label=string("BSPF, ",round(timeBSPF;digits=6),"s"))
tYLimsHere = ylims(ppp)
plot!(ones(2,1)*llKF,[tYLimsHere[1],tYLimsHere[2]],line=(:red,0.9,2,:line),label=string("KF,     ",round(timeKF;digits=6),"s"))
title!(string("LL at θ0: KF vs BSPF (M=",Mbspf," particles)"))

display(ppp)


sFileForOutput         = string("plot_mL1_BSPFcheck_M",Mbspf,".png")
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)
