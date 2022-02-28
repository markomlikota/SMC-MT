# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Create plots based on output from "main_MultipleProposalsAndRuns.jl" (so for univariate case)



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


Targetspec      = "A0"  # A0 = standard normal
SMCspec         = "1"

setOfProposals  = 1     # integer, which folder to analyze?

vIndDM          = 3   # which distance measures (DM) to consider? (see fGetDistanceMeasures.jl)



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


# cd()
# sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareSimpleExample/")
sMyPath            = pwd()



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
sSpecPath               = sMyPath * "/" * "Output/" * sSpecFolder

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

vMeasureLabsShort   = ["Area under product of pdfs","Kullback-Leibler Distance","\$\\mathcal{D}(M_0,M_1)\$","V-M_0 of pdf ratio"]
vMeasureLabs        = string.(vMeasureLabsShort," (target=N(0,1))")



# Load parameters:

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



# Load SMC-Statistics:

vStatsNames         = ["mean","std","skewness","kurtosis","p2.5","p5","p95","p97.5","'MDD'","nStages","time"]

aStdsOfStats        = zeros(nMeans0,nStds0,9)
mMeanNStages        = zeros(nMeans0,nStds0)
mMeanTimes          = zeros(nMeans0,nStds0)
mvTSs               = Array{Vector{Float64}}(undef,nMeans0,nStds0)

aDistanceMeasures   = zeros(nMeans0,nStds0,4)
[aDistanceMeasures[:,:,dd] = reshape(mDistanceMeasures[:,dd],nMeans0,nStds0) for dd=1:4]


for mm = 1:nMeans0, ss = 1:nStds0


    sProposalFileName     = string(sOutputFilePrefix,"p_m",mm,"_s",ss)


    # Precision, Runtime & #Stages:

    sReadName             = sProposalFileName * "_statstics.csv"
    # dfStatsHere           = CSV.read(sOutputPath * sReadName, DataFrame;header=1)
    # mStatsHere            = Matrix{Float64}(dfStatsHere)
    mStatsHere            = readdlm(sOutputPath * sReadName, ',', Float64, '\n'; skipstart=1)

    aStdsOfStats[mm,ss,:] = std(mStatsHere[:,1:end-2],dims=1)

    mMeanNStages[mm,ss]   = mean(mStatsHere[:,end-1])
    mMeanTimes[mm,ss]     = mean(mStatsHere[:,end])


    # Tempering Schedule:

    sReadName             = sProposalFileName * "_temperingschedule.csv"
    f                     = open(sOutputPath * sReadName,"r")
    msTSHere              = readdlm(f, ',', skipstart=1)
    replace!(msTSHere, "" => NaN)
    # NϕLongest             = size(msTSHere,2)
    # vMedianTS             = zeros(NϕLongest)
    # [vMedianTS[nn]        = mean(filter(!isnan,msTSHere[:,nn])) for nn = 1:NϕLongest]
    vMedianTS             = msTSHere[5,:]

    mvTSs[mm,ss]          = vMedianTS


end



# -- Plot distance measures vs mean/std of proposal --

vStdIndicesToUse  = [1,5,7,10]

for measureInd = vIndDM

  fGetPlotOneDistanceMeasure(mParamInfo,mDistanceMeasures,vStdIndicesToUse,measureInd)
  savefig(string("plot_distancemeasure",measureInd,"_lines.png"))

  mDistanceMeasureHere = reshape(mDistanceMeasures[:,measureInd],(nMeans0,nStds0))
  heatmap(vStds0,vMeans0,mDistanceMeasureHere,color=:reds,ylabel="μ0",xlabel="σ0",title=vMeasureLabs[measureInd])
  savefig(string("plot_distancemeasure",measureInd,"_heatmap.png"))

end

fGetPlotTwoDistanceMeasures(mParamInfo,mDistanceMeasures,vStdIndicesToUse,[2,3])
savefig(string("plot_distancemeasures23.png"))

fGetPlotTwoDistanceMeasures(mParamInfo,mDistanceMeasures,vStdIndicesToUse,[2,4])
savefig(string("plot_distancemeasures24.png"))



# -- Plot Runtime --


sTitleHere = ""#"Average Runtime [s] (target=N(0,1))"


# vs (μ0,σ0):

# plot(vStds0,vMeans0,mMeanTimes,st=:contour,color=:reds,colorbar = false)?
heatmap(vStds0,vMeans0,mMeanTimes,color=:reds,ylabel="μ0",xlabel="σ0",title=sTitleHere)
savefig("plot_meantimes_params")


# vs distance metrics:

vMeanTimes          = reshape(mMeanTimes,nMeans0*nStds0)

for measureInd = vIndDM
    myBlueHere        = cgrad(:blues)[1/nStds0]
    ppp= scatter(aDistanceMeasures[:,1,measureInd],mMeanTimes[:,1],color=myBlueHere,markerstrokewidth=0.5,markersize=10,xlabel=vMeasureLabsShort[measureInd],title=sTitleHere,label=string("σ = ",vStds0[1]),xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,legend = :topleft, background_color_legend = nothing)
    for ss = 2:nStds0
        myBlueHere        = cgrad(:blues)[ss/nStds0]
        ss < nStds0 ? sLabelHere= "" : sLabelHere = string("σ = ",vStds0[ss])
        scatter!(aDistanceMeasures[:,ss,measureInd],mMeanTimes[:,ss],color=myBlueHere,markerstrokewidth=0.5,markersize=10,label=sLabelHere)
    end
    yLimsNow = ylims(ppp)
    ylims!((0,yLimsNow[2]))
    savefig(string("plot_meantimes_distancemeasure",measureInd,"_colored"))
end



# -- #Stages --

sTitleHere = ""#"Average Nϕ (target=N(0,1))"


# vs (μ0,σ0):

heatmap(vStds0,vMeans0,mMeanNStages,color=:reds,ylabel="μ0",xlabel="σ0",title=sTitleHere)
savefig("plot_nstages_params")


# vs distance metrics:

for measureInd = vIndDM
    myBlueHere        = cgrad(:blues)[1/nStds0]
    ppp= scatter(aDistanceMeasures[:,1,measureInd],mMeanNStages[:,1],color=myBlueHere,markerstrokewidth=0.5,markersize=10,xlabel=vMeasureLabsShort[measureInd],title=sTitleHere,label=string("σ = ",vStds0[1]),xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,legend = :topleft, background_color_legend = nothing)
    for ss = 2:nStds0
        myBlueHere        = cgrad(:blues)[ss/nStds0]
        ss < nStds0 ? sLabelHere= "" : sLabelHere = string("σ = ",vStds0[ss])
        scatter!(aDistanceMeasures[:,ss,measureInd],mMeanNStages[:,ss],color=myBlueHere,markerstrokewidth=0.5,markersize=10,label=sLabelHere)
    end
    yLimsNow = ylims(ppp)
    ylims!((0,yLimsNow[2]))
    savefig(string("plot_nstages_distancemeasure",measureInd,"_colored"))
end



# -- Precision --

sTitleHere = ""#"Std of Mean (target=N(0,1))"


# vs (μ0,σ0):

heatmap(vStds0,vMeans0,aStdsOfStats[:,:,1],color=:reds,ylabel="μ0",xlabel="σ0",title=sTitleHere)
savefig("plot_std_mean_params")


# vs distance metrics:

# vStdMean          = reshape(aStdsOfStats[:,:,1],nMeans0*nStds0)
#
# for measureInd = vIndDM
#     scatter(mDistanceMeasures[:,measureInd],vStdMean,color=myBlue,markerstrokewidth=0.5,markersize=10,xlabel=vMeasureLabsShort[measureInd],title=sTitleHere,label="")
#     savefig(string("plot_std_mean_distancemeasure",measureInd))
# end

for measureInd = vIndDM
    myBlueHere        = cgrad(:blues)[1/nStds0]
    ppp= scatter(aDistanceMeasures[:,1,measureInd],aStdsOfStats[:,1,1],color=myBlueHere,markerstrokewidth=0.5,markersize=10,xlabel=vMeasureLabsShort[measureInd],title=sTitleHere,label=string("σ = ",vStds0[1]),xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,legend = :topleft, background_color_legend = nothing)
    for ss = 2:nStds0
        myBlueHere        = cgrad(:blues)[ss/nStds0]
        ss < nStds0 ? sLabelHere= "" : sLabelHere = string("σ = ",vStds0[ss])
        scatter!(aDistanceMeasures[:,ss,measureInd],aStdsOfStats[:,ss,1],color=myBlueHere,markerstrokewidth=0.5,markersize=10,label=sLabelHere)
    end
    yLimsNow = ylims(ppp)
    ylims!((0,round(yLimsNow[2];digits=1)))
    savefig(string("plot_std_mean_distancemeasure",measureInd,"_colored"))
end



# -- Plot two proposals with same distance measure but different runtimes/NStages --

myRed         = cgrad(:reds)[0.8]

vValuesForDMs = [-5.0,0.4]
vvIndicesToTake = [[1,2],[2,3]]
vMeasureLabsShortest = ["KL","areamin"]

vIndDMs = [2,3]

for mm = 1:length(vIndDMs)
    measureInd = vIndDMs[mm]
    valueHere = vValuesForDMs[mm]
    sort(abs.(mDistanceMeasures[:,measureInd] .- valueHere))
    vIndProps = sortperm(abs.(mDistanceMeasures[:,measureInd] .- valueHere))[vvIndicesToTake[mm]]

    mParamInfoHere = mParamInfo[vIndProps,:]
    vDMHere = mDistanceMeasures[vIndProps,measureInd]
    vTimesHere = vMeanTimes[vIndProps]


    dProposal1 = fGetProposalObjects(mParamInfoHere[1,3:4]...)[3]
    dProposal2 = fGetProposalObjects(mParamInfoHere[2,3:4]...)[3]

    vx            = -5:0.01:5

    plot(vx,pdf.(dTarget,vx), line =(:black,0.9,2, :line), label="target")
    plot!(vx,pdf.(dProposal1,vx), line =(myRed,0.9,2, :line), label=string(vMeasureLabsShortest[mm],"=",round(vDMHere[1];digits=4),", time=",round(vTimesHere[1];digits=4)))
    plot!(vx,pdf.(dProposal2,vx), line =(myBlue,0.9,2, :line), label=string(vMeasureLabsShortest[mm],"=",round(vDMHere[2];digits=4),", time=",round(vTimesHere[2];digits=4)))

    savefig(string("plot_twoprops_distancemeasure",measureInd))
end



# -- Tempering Schedules --

sTitleHere      = "Median Tempering Schedules"

vMeansToConsider = vMeans0[1:1:nMeans0]
nMeansToConsider = length(vMeansToConsider)

for ss = 1:nStds0
    myBlueHere        = cgrad(:blues)[1/nMeansToConsider]
    vTShere           = mvTSs[1,ss]
    ppp= plot(1:length(vTShere),vTShere,line =(myBlueHere,0.9,2, :line),xlabel="stage n",title="",label=string("μ=",vMeansToConsider[1]),xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,legend = :bottomright, background_color_legend = nothing)
    for mm = 2:nMeansToConsider
        myBlueHere        = cgrad(:blues)[mm/nMeansToConsider]
        vTShere           = mvTSs[mm,ss]
        mm < nMeansToConsider ? sLabelHere= "" : sLabelHere = string("μ=",vMeansToConsider[mm])
        plot!(1:length(vTShere),vTShere,line =(myBlueHere,0.9,2, :line),label=sLabelHere)
    end
    savefig(string("plot_TS_ss",ss,".png"))
end


vStdsToConsider = vStds0[1:1:nStds0]
nStdsToConsider = length(vStdsToConsider)

for mm = 1:nMeans0
    myBlueHere        = cgrad(:blues)[1/nStdsToConsider]
    vTShere           = mvTSs[mm,1]
    ppp= plot(1:length(vTShere),vTShere,line =(myBlueHere,0.9,2, :line),xlabel="stage n",title="",label=string("σ=",vStds0[1]),xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,legend = :bottomright, background_color_legend = nothing)
    for ss = 2:nStdsToConsider
        myBlueHere        = cgrad(:blues)[ss/nStdsToConsider]
        vTShere           = mvTSs[mm,ss]
        ss < nStdsToConsider ? sLabelHere= "" : sLabelHere = string("σ=",vStds0[ss])
        plot!(1:length(vTShere),vTShere,line =(myBlueHere,0.9,2, :line),label=sLabelHere)
    end
    savefig(string("plot_TS_mm",mm,".png"))
end

#
# # Mock plot to get theta printed for figure 1 of paper:
#
# plot(collect(1:10),ones(10),xlabel="\$\\theta\$")
# savefig(string("plot_theta.png"))
