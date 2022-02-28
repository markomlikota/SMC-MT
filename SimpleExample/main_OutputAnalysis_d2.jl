# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Create plots based on output from "main_MultipleProposalsAndRuns_d2.jl" (so for bivariate case)



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


Targetspec      = "B0"  # A0 = standard normal
SMCspec         = "1"

setOfProposals  = 2     # integer, which folder to analyze?

vIndDM          = 1:4   # which distance measures (DM) to consider? (see fGetDistanceMeasures.jl)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------



using Plots
pyplot()
using Roots
using StatsBase
using KernelDensity
using DataFrames
using CSV
using StatsFuns
using DelimitedFiles
using Distributions
using LinearAlgebra



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


cd()
sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareSimpleExample/")



include(sMyPath * "/Functions/fHelpers.jl")
include(sMyPath * "/Functions/fFolderFileManagement.jl")
include(sMyPath * "/Functions/fMyPlot.jl")
include(sMyPath * "/Functions/fGetProposalObjects.jl")
include(sMyPath * "/Functions/fGetDistanceMeasures.jl")
include(sMyPath * "/Functions/fFunctionsOutputAnalysis.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------



# -- Preliminary ---------------------------------------------------------------


# Set up directory for estimation output:

sSpecFolder             = string("t",Targetspec,"_smc",SMCspec,"/")
sSpecPath               = sMyPath * "Output/" * sSpecFolder

sHelperFolder           = string("SetOfProposals",setOfProposals) * "/"
sHelperPath             = sSpecPath * sHelperFolder

sEstSpecFolder          = ""
sEstSpecPath            = sHelperPath * sEstSpecFolder


# Set up target density and SMC settings:

include(sMyPath * "/SpecFiles/" * string("script_Targetspec",Targetspec,".jl") )

include(sMyPath * "/SpecFiles/" * string("script_SMCspec",SMCspec,".jl") )


# Set up output folder and filename prefix:

sOutputPath       = sEstSpecPath
sOutputFilePrefix = string("t",Targetspec,"_smc",SMCspec,"_")



# -- Plots ---------------------------------------------------------------------


cd(sOutputPath)

myBlue              = cgrad(:blues)[1.0]

vMeasureLabsShort   = ["Area under product of pdfs","Kullback-Leibler Distance","Area under minimum of pdfs","V-M_0 of pdf ratio"]
vMeasureLabs        = string.(vMeasureLabsShort," (target= 2 x N(0,1), ρ=0.25)")



# Load parameters:

sReadName            = sOutputFilePrefix * "parameters.csv"
# dfParams             = CSV.read(sOutputPath * sReadName, DataFrame;header=1)
# mParams              = Matrix{Float64}(dfParams)
mParams              = readdlm(sOutputPath * sReadName, ',', Float64, '\n'; skipstart=1)

mParamInfo           = mParams[:,1:10]
mDistanceMeasures    = mParams[:,11:end]
nProps               = size(mParamInfo,1)

vMeans01             = unique(mParamInfo[:,6])
vStds01              = unique(mParamInfo[:,7])
vMeans02             = unique(mParamInfo[:,8])
vStds02              = unique(mParamInfo[:,9])
vRhos0               = unique(mParamInfo[:,10])

nMeans01             = length(vMeans01)
nStds01              = length(vStds01)
nMeans02             = length(vMeans02)
nStds02              = length(vStds02)
nRhos0               = length(vRhos0)



# Load SMC-Statistics:

vStatsNames         = ["mean","std","skewness","kurtosis","p2.5","p5","p95","p97.5","'MDD'","nStages","time"]

mStdsOfStats        = zeros(nProps,9)
vMeanNStages        = zeros(nProps)
vMeanTimes          = zeros(nProps)
vvTSs               = Array{Vector{Float64}}(undef,nProps)


for kk = 1:nProps


    mm1,ss1,mm2,ss2,pp    = Int.(mParamInfo[kk,1:5])
    μ1g,σ1g,μ2g,σ2g,ρg    = mParamInfo[kk,6:10]


    sProposalFileName     = sOutputFilePrefix * "p_" * fCombinevStrings(Int.(mParamInfo[kk,1:5]))


    # Precision, Runtime & #Stages:

    sReadName             = sProposalFileName * "_statstics.csv"
    # dfStatsHere           = CSV.read(sOutputPath * sReadName, DataFrame;header=1)
    # mStatsHere            = Matrix{Float64}(dfStatsHere)
    mStatsHere            = readdlm(sOutputPath * sReadName, ',', Float64, '\n'; skipstart=1)

    mStdsOfStats[kk,:]    = std(mStatsHere[:,1:end-2],dims=1)

    vMeanNStages[kk]      = mean(mStatsHere[:,end-1])
    vMeanTimes[kk]        = mean(mStatsHere[:,end])


    # Tempering Schedule:

    sReadName             = sProposalFileName * "_temperingschedule.csv"
    f                     = open(sOutputPath * sReadName,"r")
    msTSHere              = readdlm(f, ',', skipstart=1)
    replace!(msTSHere, "" => NaN)
    NϕLongest             = size(msTSHere,2)
    vMedianTS             = zeros(NϕLongest)
    [vMedianTS[nn]        = median(filter(!isnan,msTSHere[:,nn])) for nn = 1:NϕLongest]

    vvTSs[kk]             = vMedianTS


end




# -- Plot Runtime --

sTitleHere = "Average Runtime [s] (target=N(0,1))"


# vs distance metrics:

for measureInd = vIndDM
    scatter(mDistanceMeasures[:,measureInd],vMeanTimes,color=myBlue,markerstrokewidth=0.5,markersize=10,xlabel=vMeasureLabsShort[measureInd],title=sTitleHere,label="")
    savefig(string("plot_meantimes_distancemeasure",measureInd))
end



# -- #Stages --

sTitleHere = "Average Nϕ (target=N(0,1))"


# vs distance metrics:

for measureInd = vIndDM
    scatter(mDistanceMeasures[:,measureInd],vMeanNStages,color=myBlue,markerstrokewidth=0.5,markersize=10,xlabel=vMeasureLabsShort[measureInd],title=sTitleHere,label="")
    savefig(string("plot_nstages_distancemeasure",measureInd))
end



# -- Precision --

sTitleHere = "Std of Mean (target=N(0,1))"


# vs distance metrics:

vStdMean          = mStdsOfStats[:,1]

for measureInd = vIndDM
    scatter(mDistanceMeasures[:,measureInd],vStdMean,color=myBlue,markerstrokewidth=0.5,markersize=10,xlabel=vMeasureLabsShort[measureInd],title=sTitleHere,label="")
    savefig(string("plot_std_mean_distancemeasure",measureInd))
end
