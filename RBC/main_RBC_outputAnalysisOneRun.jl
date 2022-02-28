# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Creates plots of estimation output of one estimation specification (DGP-Model-Prior-SMC; for LT or MT) and one run.

# - (marginal) posterior and stage1-density (prior/M0-posterior(w./w.o. early stopping))
# - (marginal) posterior



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec         = "3"     # data from which DGP to take?
Model1spec      = "L2"    # which model (LL) to estimate?
Prior1spec      = "L22"   # which prior to take?
SMC1spec        = "03"    # which values for SMC tuning parameters to take?

sTemperingType  = "LT"          # "LT" or "MT"; how to estimate M1?
ϕ_Nϕ            = 1.0           # stop early?


# Only relevant if sTemperingType == "MT":

Model0spec      = "L1"            # which model is M0?
Prior0spec      = "L12"           # which prior was taken for M0?
ϕ_Nϕ_M0         = 0.5             # early stopping for M0?
SMC0spec        = "03"    # which values for SMC tuning parameters were taken for M0?



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


using Plots
pyplot()
using StatsBase
using KernelDensity
using DataFrames
using CSV
using DelimitedFiles
using LinearAlgebra
using Statistics



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


cd()
sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/RBC/")
# sMyPath            = pwd() * "/"


include(sMyPath * "Functions/fHelpersSMC.jl")
include(sMyPath * "Functions/fFolderFileManagement.jl")
include(sMyPath * "Functions/fFunctionsOutputAnalysis.jl")
include(sMyPath * "Functions/fPrPoPlots.jl")
include(sMyPath * "Functions/fMyPlot.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# -- Preliminary ---------------------------------------------------------------


sSpecPath       = sMyPath * "Output/" * fGetSpecFolder(DGPspec,Model1spec,Prior1spec,SMC1spec)

if sTemperingType == "LT"
    sFilePath       = sSpecPath * fGetEstSpecFolderLT(ϕ_Nϕ)
else
    sFilePath       = sSpecPath * fGetEstSpecFolderMT(ϕ_Nϕ_M0,Model0spec,Prior0spec,SMC0spec)
end
sPlotPath       = sFilePath

myBlue          = cgrad(:blues)[1.0]
myRed           = cgrad(:reds)[0.8]

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )


# -- Posterior Evolution & Prior-Posterior Plots -------------------------------


# Object with Posterior Particles:

sFileSuffix             = "StageLast_Particles.csv"

mParticlesStageLast     = fGetFile(sFilePath,sFileSuffix)


# Object with Prior Particles:

sFileSuffix             = "Stage1_Particles.csv"

mParticlesStage1        = fGetFile(sFilePath,sFileSuffix)


# Combine into one object:

N,nP                    = size(mParticlesStage1)

aAllParticles           = Array{Float64}(undef,N,nP,2)
aAllParticles[:,:,1]    = mParticlesStage1
aAllParticles[:,:,2]    = mParticlesStageLast


# Plots:

vBoundTypes = [2,1,2,1,1,1,2,2,1,1,0,1,1,1] # for PrPoPlots; based on parameter domains defined in fIsDrawValid in Priorspec

vDoesParamVary  = [length(unique(aAllParticles[:,pp,1]))!=1 for pp = 1:nP]

vPIndVary       = collect(1:nP)[vDoesParamVary]


for pInd in vPIndVary

    # Preliminary:

    yLab        = vParLabs[pInd]
    mWeights    = ones(N,9)
    boundsType  = Int(vBoundTypes[pInd])


    # Plot posterior and stage1-density (prior or M0-posterior(w./w.o. early stopping)):

    ppp = fPrPoPlot(pInd,aAllParticles,mWeights,boundsType,myBlue)
    title!(vParLabs[pInd])

    vzLimsHere = ylims(ppp)

    # Add red vertical line at true value:
    plot!(vec(ones(2,1))*θ0[pInd], [vzLimsHere[1],vzLimsHere[2]],line=(myRed,3,:dot),label="")

    savefig(sPlotPath * string("plot_prpo_p",pInd,".png"))


    # Plot only posterior:

    kPo         = kde(aAllParticles[:,pInd,end], weights = Weights(mWeights[:,end]))
    stPo, enPo  = fFindPlotBounds(boundsType,kPo.x)
    ppp         = fMyPlot(kPo.x[stPo:enPo],kPo.density[stPo:enPo],vParLabs[pInd],"","","",myBlue)
    title!(vParLabs[pInd])

    vzLimsHere = ylims(ppp)

    # Add red vertical line at true value:
    plot!(vec(ones(2,1))*θ0[pInd], [vzLimsHere[1],vzLimsHere[2]],line=(myRed,3,:dot),label="",yticks=false)

    savefig(sPlotPath * string("plot_po_p",pInd,".png"))

end



# -- Prior HPD -----------------------------------------------------------------

using StatisticalRethinking

vvHPDbounds         = [hpdi(mParticlesStage1[:,pp]) for pp = 1:nP]

mHPDbounds          = zeros(nP,2)
[mHPDbounds[pp,:]   = vvHPDbounds[pp] for pp = 1:nP]

# # Note that it's not the same as 5th and 95th percentile (last two columns here):
# using Distributed
# include(sMyPath * "Functions/fSetupSMC.jl")
# mPriorStats         = fGetParticleStats(ones(N),mParticlesStage1)
# [mPriorStats mHPDbounds]

[vParLabs mHPDbounds]
