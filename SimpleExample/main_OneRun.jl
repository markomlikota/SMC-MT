# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Runs SMC one time given target and proposal distribution (& SMC settings).



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


Targetspec  = "A0"
Propspec    = "A0"
SMCspec     = "1"



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



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------



# -- Set up directory for estimation output ------------------------------------


sSpecFolder             = string("t",Targetspec,"_p",Propspec,"_smc",SMCspec,"/")
sHelperFolder           = ""
sEstSpecFolder          = ""


sSpecPath       = sMyPath * "Output/" * sSpecFolder
sHelperPath     = sSpecPath * sHelperFolder
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


# Set up proposal and target densities and SMC settings:

@everywhere @eval Propspec    = $Propspec
@everywhere @eval Targetspec  = $Targetspec

@everywhere include(sMyPath * "/SpecFiles/" * string("script_Propspec",Propspec,".jl") )
@everywhere include(sMyPath * "/SpecFiles/" * string("script_Targetspec",Targetspec,".jl") )

include(sMyPath * "/SpecFiles/" * string("script_SMCspec",SMCspec,".jl") )


# Setup output folder and filename prefix:

sOutputPath       = sEstSpecPath
sOutputFilePrefix = string("t",Targetspec,"_p",Propspec,"_smc",SMCspec,"_")


# Define functions/objects for SMC:

@everywhere begin

  fLLtilde(θ)         = fEvalLogProposal(θ)
  fLL(θ)              = fEvalLogTarget(θ)
  #fPriorDraw()        = fPriorDrawM1() #not needed here; (all) initial particles supplied
  fIsDrawValid(θ)     = 1
  fPriorLogEval(θ)    = 0

end

vParLabs  = string.("p",1:nP)
@everywhere @eval vParLabs    = $vParLabs


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

setSeed               = 124

vmParticles, vvWeights, vmLogLiks, vvPriors, vAcceptRate, vSstep, vESS, vc, vϕ, vLogMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed,mInitParticles,vPIndDraw)



# -- Visualization -------------------------------------------------------------


cd(sOutputPath)


# Preliminary:

# Delete all plots in directory: (number of stages in new run could be lower than in previous, so taht not all plots get overwritten)
foreach(rm, filter(endswith(".png"), readdir(sOutputPath,join=true)))


nPlotStages   = length(vmParticles)

myBlue        = cgrad(:blues)[0.8]
myBlue2       = cgrad(:blues)[0.2]

myRed         = cgrad(:reds)[0.8]
myRed2        = cgrad(:reds)[0.2]

vDarkReds     = cgrad([cgrad(:reds)[0.2],cgrad(:reds)[1.0]])
vDarkBlues    = cgrad([cgrad(:blues)[0.2],cgrad(:blues)[1.0]])
vColors       = range(myRed2, stop=myBlue2, length=nPlotStages)

if nP == 1 # define nice x-axis (to be held constant) in case of univariate dist.
    vx            = round(minimum([mInitParticles[:,1]...,-5]);digits=1):0.1:round(maximum([mInitParticles[:,1]...,5]);digits=1)
else
    vx            = round(minimum([vec(mInitParticles)...,-5]);digits=1):0.1:round(maximum([vec(mInitParticles)...,5]);digits=1)
end

vStagesToPlot = Int.(round.(range(1,nPlotStages,length=nPlotStages)))


# Plots:

if nP == 1

  for nn in vStagesToPlot

    histCol = vColors[nn]

    histogram(vmParticles[nn][:,1],label=string("particles"),weights=vvWeights[nn],linecolor=histCol,normalize=:pdf,fillcolor=histCol,background_color_legend = nothing,xtickfont=font(14),ytickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14)

    plot!(vx,pdf.(dProposal,vx), line =(myRed,0.9,2, :line), label="proposal")
    plot!(vx,pdf.(dTarget,vx), line =(myBlue,0.9,2, :line), label="target")

    #title!(string("Distribution Evolution; Stage ",nn-1,", ϕ =",round(vϕ[nn];digits=4)))

    savefig(string("plot_stage",nn-1,".png"))

  end

elseif nP == 2

  for nn in vStagesToPlot

    histCol = vColors[nn]

    scatter(vmParticles[nn][:,2],vmParticles[nn][:,1],marker = (histCol),markerstrokewidth=0.1,markersize=vvWeights[nn]*2,legend=:none)

    mPDFg = [pdf(dProposal,[i,j]) for i in vx, j in vx]
    plot!(vx,vx,mPDFg,st=:contour,color=vDarkReds,linewidth=2,label="proposal",colorbar = false,levels=5)

    mPDFf = [pdf(dTarget,[i,j]) for i in vx, j in vx]
    plot!(vx,vx,mPDFf,st=:contour,color=vDarkBlues,linewidth=2,label="target",colorbar = false,levels=5)

    title!(string("Distribution Evolution; Stage ",nn-1,", ϕ =",round(vϕ[nn];digits=4)))

    savefig(string("plot_stage",nn-1,".png"))

  end

end
