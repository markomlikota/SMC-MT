# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Script to check a few things for model "VFI"



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



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


@everywhere cd()
@everywhere sMyPath            = string(pwd(),"/Dropbox/ModelTempering/Software/RBC/")
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

@everywhere include(sMyPath * "SpecFiles/" * string("script_ModelspecVFI.jl") )


# to get θ0:

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )



# -- Compare Simulated Data under mL1 and mVFI ---------------------------------


sFolderForOutput    = sMyPath * "Data/"

mData_mL1   = fGetFile(sMyPath * "Data/","Data_" * string("dgp",DGPspec) * "_mL1" * ".csv")
mData_mVFI  = fGetFile(sMyPath * "Data/","Data_" * string("dgp",DGPspec) * "_mVFI" * ".csv")

myBlue          = cgrad(:blues)[1.0]
myRed           = cgrad(:reds)[0.8]

vYlabs      = ["log y^o","log i^o","log l^o"]

p1          = fMyPlot(1:T,mData_mL1[:,1],vYlabs[1],"","","L1",myBlue)
plot!(1:T,mData_mVFI[:,1],label="VFI",line=(myRed,0.9,2,:line))

p2          = fMyPlot(1:T,mData_mL1[:,2],vYlabs[2],"","","L1",myBlue)
plot!(1:T,mData_mVFI[:,2],label="VFI",line=(myRed,0.9,2,:line))

p3          = fMyPlot(1:T,mData_mL1[:,3],vYlabs[3],"t","","L1",myBlue)
plot!(1:T,mData_mVFI[:,3],label="VFI",line=(myRed,0.9,2,:line))

pp          = plot(p1,p2,p3,layout=(3,1))
display(pp)

sFileForOutput         = "Data_" * string("dgp",DGPspec) *"_L1vsVFI" * ".png"
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)



# -- Investigate conditional distribution of investment given past states ------


sFolderForOutput    = sMyPath * "Output/Checks/"

θ       = deepcopy(θ0)

fSol_K_, maxDifference, iteration = fSolveModelVFI(θ,1)

ρz, σz, ρb, σb, r, τ, α, δ, ν, ϕ1, ϕ2, σy, σi, σl = θ

cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ)
#
# fI(fSol_K_(kSS,Zstar,Bstar),kSS,δ)

# _vs = [kSS, lSS, ySS, iSS]
# _vε = [0,0]
# vε = fFε(_vε,ρz,σz,ρb,σb)


nk      = 50
nMC     = 1#100000
# Note that mean of I and log I below change every time this code is run (esp. for low nMC). One could actually also just take the means to be Ihere and logIhere and add the bands around them...

mInvestment     = zeros(nMC,nk)
mLogInvestment  = zeros(nMC,nk)

fracEitherSide = 0.025
vkgrid  = fGetKgrid(kSS,nk,fracEitherSide)

vε = [0,0]

for kk = 1:nk

    _vsHere = [vkgrid[kk],lSS,ySS,iSS]

    for ee = 1:nMC

        Ihere       = fΦ(_vsHere,vε,fSol_K_,θ,Zstar,Bstar)[4]
        logIhere    = log(Ihere)

        vη          = 0 #rand(Normal(0,σi))

        mInvestment[ee,kk]      = Ihere + vη
        mLogInvestment[ee,kk]   = logIhere + vη

    end

end

vMeanI  = mean(mInvestment,dims=1)'
vStdI   = σi * ones(nk)

pp = fMyRibbonPlot(vkgrid,vMeanI,1.64*vStdI,"p(I|K,Z*,B*)","K","")
display(pp)

sFileForOutput         = string("plot_mVFI_dist_I_dgp",DGPspec,".png")
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)


# Log I:

vMeanLogI  = mean(mLogInvestment,dims=1)'
# vStdI   = σi * ones(nk)

pp = fMyRibbonPlot(vkgrid,vMeanLogI,1.64*vStdI,"p(log I|K,Z*,B*)","K","")
display(pp)

sFileForOutput         = string("plot_mVFI_dist_logI_dgp",DGPspec,".png")
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)



# -- Investigate whether VFI converges for many prior draws --------------------


sFolderForOutput    = sMyPath * "Output/Checks/"


# This is now done inside "Priorspec"-functions, but keep it here just in case...:
# @everywhere @eval DGPspec        = $DGPspec # we need DGPspec, because prior fixes some parameters at true values, which are contained in DGPspec
# @everywhere include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )

@everywhere @eval DGPspec        = $DGPspec # we need DGPspec, because prior fixes some parameters at true values, which are contained in DGPspec


@everywhere include(sMyPath * "SpecFiles/" * string("script_PriorspecVFI1.jl") )


nDraws = 1000
vMaxDiffs   = SharedArray{Float64}((nDraws))
mPriorDraws = SharedArray{Float64}((nDraws,nP))
@sync @distributed for dd = 1:nDraws

    println(dd)

    θ       = fPriorDraw_VFI1()

    mPriorDraws[dd,:] = θ

    fSol_K_, maxDiff, iteration = fSolveModelVFI(θ,1)

    vMaxDiffs[dd] = maxDiff

end
