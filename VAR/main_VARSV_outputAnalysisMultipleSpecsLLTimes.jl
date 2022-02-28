# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Creates plots of estimation output of several DGP-SMC-specifications for all estimation specifications (M1-LT, M1-MT) over many runs.



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


vModel1specs            = ["SV","SVx2","SVx4","SVx8"]
vPrior1specs            = ["SV1","SVx21","SVx41","SVx81"]

DGPspec                 = "1"     # data from which DGP to take?
SMCspec                 = "1"     # which values for SMC tuning parameters to take?

Model0spec              = "HM"            # which model is M0?
Prior0spec              = "HM1"           # which prior to take for M0?

vPhiLast                = [0.2, 0.4, 0.6, 0.8, 1.0]   # (0.0 automatically considered)

# Remaining options (see main_VARSV_OneRun.jl) specified inside functions to extract files in fFunctionsOutputAnalysis.jl

nRun                = 1

sRunFolderName      = "xrun"



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


# cd()
# sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareVAR/")
sMyPath             = pwd() * "/"


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


sPlotPath       = sMyPath * "Output/"

nPhiLast        = length(vPhiLast)

myBlue          = cgrad(:blues)[1.0]

vMyShapes       = [:solid, :dash, :dashdot, :dot]

vXaxis          = [0,vPhiLast...]



# -- Obtain Time & MDD across all specs ----------------------------------------


mSpecs          = [vModel1specs vPrior1specs]
nSpecs          = size(mSpecs,1)

mTimes          = zeros(1+nPhiLast,nSpecs)
mTimesStd       = zeros(1+nPhiLast,nSpecs)
# mTimesLower     = zeros(1+nPhiLast,nSpecs)
# mTimesUpper     = zeros(1+nPhiLast,nSpecs)
# mLogMDDStd      = zeros(1+nPhiLast,nSpecs)

for sspec in 1:nSpecs

    local Model1specHere   = mSpecs[sspec,1]
    local Prior1specHere   = mSpecs[sspec,2]
    sFileSuffix     = "FinalStats.csv"
    vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast,sRunFolderName,Model1specHere,Prior1specHere)

    vTimes          = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"mean")
    vTimesStds      = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"std")
    # vTimesLower     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"lower")
    # vTimesUpper     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"upper")
    # vLogMDDStd      = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"std")

    mTimes[:,sspec]        = vTimes
    mTimesStd[:,sspec]     = vTimesStds
    # mTimesLower[:,sspec]   = vTimesLower
    # mTimesUpper[:,sspec]   = vTimesUpper
    # mLogMDDStd[:,sspec]    = vLogMDDStd

end

mTimes          = mTimes./60
mTimesStd       = mTimesStd./60
# mTimesLower     = mTimesLower./60
# mTimesUpper     = mTimesUpper./60



# -- Time & MDD plots across all specs -----------------------------------------


vLineLabs       = vModel1specs

xLabPhiLast     = "\$\\psi_*\$" # "\$\\phi_{N_\\phi}(M_0)\$"


# Times Plot:

plotTimes       = plot(vXaxis,mTimes[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="minutes",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mTimes[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mTimes[:,3],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])
plot!(vXaxis,mTimes[:,4],line=(myBlue,0.9,2,vMyShapes[4]),label=vLineLabs[4])

xticks!(vXaxis)
ylims!((0,ylims(plotTimes)[2]))

savefig(sPlotPath * "plot_meanTimes_LLtimes.png")


# Relative Times Plot:

mRelTimes    = mTimes./mTimes[1,:]'

plotRelTimes    = plot(vXaxis,mRelTimes[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="minutes relative to LT",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mRelTimes[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mRelTimes[:,3],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])
plot!(vXaxis,mRelTimes[:,4],line=(myBlue,0.9,2,vMyShapes[4]),label=vLineLabs[4])

xticks!(vXaxis)
ylims!((0,ylims(plotRelTimes)[2]))

savefig(sPlotPath * "plot_meanRelTimes_LLtimes.png")
