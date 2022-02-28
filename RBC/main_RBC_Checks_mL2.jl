# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Jan 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Script to check a few things for model "L2"



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec         = "1"     # data from which DGP to take?



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
@everywhere sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/RBC/")
# @everywhere sMyPath            = pwd() * "/"


@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")
include(sMyPath * "/Functions/fMyPlot.jl")
include(sMyPath * "/Functions/fFolderFileManagement.jl")
@everywhere include(sMyPath * "/Functions/fHelpersSMC.jl")


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


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



# -- Compare Simulated Data under mL1 and mL2 ---------------------------------


sFolderForOutput    = sMyPath * "Data/"

mData_mL1   = fGetFile(sMyPath * "Data/","Data_" * string("dgp",DGPspec) * "_mL1" * ".csv")
mData_mL2   = fGetFile(sMyPath * "Data/","Data_" * string("dgp",DGPspec) * "_mL2" * ".csv")

myBlue          = cgrad(:blues)[1.0]
myRed           = cgrad(:reds)[0.8]

vYlabs      = ["log y^o","log i^o","log l^o"]

p1          = fMyPlot(1:T,mData_mL1[:,1],vYlabs[1],"","","L1",myBlue)
plot!(1:T,mData_mL2[:,1],label="L2",line=(myRed,0.9,2,:line))

p2          = fMyPlot(1:T,mData_mL1[:,2],vYlabs[2],"","","L1",myBlue)
plot!(1:T,mData_mL2[:,2],label="L2",line=(myRed,0.9,2,:line))

p3          = fMyPlot(1:T,mData_mL1[:,3],vYlabs[3],"t","","L1",myBlue)
plot!(1:T,mData_mL2[:,3],label="L2",line=(myRed,0.9,2,:line))

pp          = plot(p1,p2,p3,layout=(3,1))
display(pp)

sFileForOutput         = "Data_" * string("dgp",DGPspec) *"_L1vsL2" * ".png"
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)
