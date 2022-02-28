# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Creates plots of estimation output of several DGP-SMC-specifications for all estimation specifications (M1-LT, M1-MT) over many runs.



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


vDGPspecs         = ["1","2","3","1"]     # data from which DGPs to take?
vSMCspecs         = ["1","1","1","2"]     # which SMC specs to take?

vPhiLast          = [0.2, 0.4, 0.6, 0.8, 1.0]   # (0.0 automatically considered)

# Remaining options (see main_VARSV_OneRun.jl) specified inside functions to extract files in fFunctionsOutputAnalysis.jl

nRun              = 200



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

vMyShapes       = [:solid, :dash, :dot]

vXaxis          = [0,vPhiLast...]



# -- Obtain Time & MDD across all specs ----------------------------------------


mSpecs          = [vDGPspecs vSMCspecs]
nSpecs          = size(mSpecs,1)

mTimes          = zeros(1+nPhiLast,nSpecs)
mTimesStd       = zeros(1+nPhiLast,nSpecs)
mTimesLower     = zeros(1+nPhiLast,nSpecs)
mTimesUpper     = zeros(1+nPhiLast,nSpecs)
mLogMDDStd      = zeros(1+nPhiLast,nSpecs)

for sspec in 1:nSpecs

    local DGPspec   = mSpecs[sspec,1]
    local SMCspec   = mSpecs[sspec,2]
    sFileSuffix     = "FinalStats.csv"
    vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)

    vTimes          = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"mean")
    vTimesStds      = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"std")
    vTimesLower     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"lower")
    vTimesUpper     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"upper")
    vLogMDDStd      = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"std")

    mTimes[:,sspec]        = vTimes
    mTimesStd[:,sspec]     = vTimesStds
    mTimesLower[:,sspec]   = vTimesLower
    mTimesUpper[:,sspec]   = vTimesUpper
    mLogMDDStd[:,sspec]    = vLogMDDStd

end

mTimes          = mTimes./60
mTimesStd       = mTimesStd./60
mTimesLower     = mTimesLower./60
mTimesUpper     = mTimesUpper./60



# -- Time & MDD plots across all DGP specs -------------------------------------


mTimesDGP       = mTimes[:,vSMCspecs.=="1"]
mTimesStdDGP    = mTimesStd[:,vSMCspecs.=="1"]
mLogMDDStdDGP   = mLogMDDStd[:,vSMCspecs.=="1"]
mRelTimesDGP    = mTimesDGP./mTimesDGP[1,:]'

vLineLabs       = ["DGP 1", "DGP 2", "DGP 3"]

xLabPhiLast     = "\$\\psi_*\$" # "\$\\phi_{N_\\phi}(M_0)\$"


# Times Plot:

plotTimes       = plot(vXaxis,mTimesDGP[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="minutes",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mTimesDGP[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mTimesDGP[:,3],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])

xticks!(vXaxis)
ylims!((0,ylims(plotTimes)[2]))

savefig(sPlotPath * "plot_meanTimes_DGPs.png")


# Relative Times Plot:

plotRelTimes    = plot(vXaxis,mRelTimesDGP[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="minutes relative to LT",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mRelTimesDGP[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mRelTimesDGP[:,3],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])

xticks!(vXaxis)
ylims!((0,ylims(plotRelTimes)[2]))

savefig(sPlotPath * "plot_meanRelTimes_DGPs.png")


# Log MDD Std Plot:

plotStds        = plot(vXaxis,mLogMDDStdDGP[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mLogMDDStdDGP[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mLogMDDStdDGP[:,3],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])

xticks!(vXaxis)
# ylims!((0,ylims(plotStds)[2]))

savefig(sPlotPath * "plot_stdLogMDD_DGPs.png")



# -- Time & MDD plots across all SMC specs -------------------------------------


mTimesSMC       = mTimes[:,vDGPspecs.=="1"]
mTimesStdSMC    = mTimesStd[:,vDGPspecs.=="1"]
mTimesLowerSMC  = mTimesLower[:,vDGPspecs.=="1"]
mTimesUpperSMC  = mTimesUpper[:,vDGPspecs.=="1"]
mLogMDDStdSMC   = mLogMDDStd[:,vDGPspecs.=="1"]
mRelTimesSMC    = mTimesSMC./mTimesSMC[1,:]'

vLineLabs       = ["N = 500","N = 1000"]


# Times Plot:

plotTimes       = plot(vXaxis,mTimesSMC[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="minutes",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mTimesSMC[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])

xticks!(vXaxis)
ylims!((0,ylims(plotTimes)[2]))

savefig(sPlotPath * "plot_meanTimes_SMCs.png")


# Times Plot with ribbon:

plotTimesRibbon = plot(vXaxis,mTimesSMC[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="minutes",title="", ribbon=(mTimesSMC[:,1]-mTimesLowerSMC[:,1],mTimesUpperSMC[:,1]-mTimesSMC[:,1]),fillalpha=0.2,fillcolor=myBlue, xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mTimesSMC[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2], ribbon=(mTimesSMC[:,2]-mTimesLowerSMC[:,2],mTimesUpperSMC[:,2]-mTimesSMC[:,2]),fillalpha=0.2,fillcolor=myBlue)

xticks!(vXaxis)
ylims!((0,ylims(plotTimesRibbon)[2]))

savefig(sPlotPath * "plot_meanTimes_SMCs_ribbon.png")



# Relative Times Plot:

plotRelTimes    = plot(vXaxis,mRelTimesSMC[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="minutes relative to LT",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mRelTimesSMC[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])

xticks!(vXaxis)
ylims!((0,ylims(plotRelTimes)[2]))

savefig(sPlotPath * "plot_meanRelTimes_SMCs.png")


# Log MDD Std Plot:

plotStds        = plot(vXaxis,mLogMDDStdSMC[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mLogMDDStdSMC[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])

xticks!(vXaxis)
# ylims!((0,ylims(plotStds)[2]))

savefig(sPlotPath * "plot_stdLogMDD_SMCs.png")



# -- std log MDD plot across all specs (DGP & SMC) -----------------------------


vLineLabs       = ["DGP 1, N = 500","DGP 2, N = 500","DGP 3, N = 500","DGP 1, N = 1000"]

vMyShapes       = [:dashdot, :solid, :dash, :dot]

# Log MDD Std Plot:

plotStds        = plot(vXaxis,mLogMDDStd[:,1],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mLogMDDStd[:,2],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mLogMDDStd[:,3],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])
plot!(vXaxis,mLogMDDStd[:,4],line=(myBlue,0.9,2,vMyShapes[4]),label=vLineLabs[4])

xticks!(vXaxis)
# ylims!((0,ylims(plotStds)[2]))

savefig(sPlotPath * "plot_stdLogMDD.png")



# -- Data (SV) plot across all specs (DGP) -------------------------------------


mdDGP1      = fGetFile(sMyPath * "Data/","SV_dgp1.csv")
mdDGP2      = fGetFile(sMyPath * "Data/","SV_dgp2.csv")
mdDGP3      = fGetFile(sMyPath * "Data/","SV_dgp3.csv")

vmd         = Matrix{Float64}[mdDGP1, mdDGP2, mdDGP3]

T           = size(vmd[1],1) -1

p1          = fMyPlot(0:T,vmd[3][:,1],"d1","")
yLimD1      = ylims(p1)
p2          = fMyPlot(0:T,vmd[3][:,2],"d2","t")
yLimD2      = ylims(p2)


for ddgp = 1:3

    p1          = fMyPlot(0:T,vmd[ddgp][:,1],"d1","")
    ylims!((0,yLimD1[2]))
    p2          = fMyPlot(0:T,vmd[ddgp][:,2],"d2","t")
    ylims!((0,yLimD2[2]))

    pp          = plot(p1,p2,layout=(2,1),legend=false)

    savefig(sMyPath * "Data/SV_dgp"*string(ddgp)*"_sameaxes.png")

end



# -- Obtain Importance Sampling Weights for all estimation specifications ------


vLineLabs       = ["DGP 1", "DGP 2", "DGP 3"]

vMyShapes       = [:solid, :dash, :dot]

sFileSuffix     = "Stage1_WLLP.csv"

aVarISweights   = zeros(nRun,nPhiLast+1,3)

for ddgp = 1:3

    DGPspec = vDGPspecs[ddgp]
    SMCspec = vSMCspecs[ddgp]

    vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)


    mVarISweights   = zeros(nRun,nPhiLast+1)

    for rr = 1:nRun


        vLL                 = vvvaStats[3][1][rr][:,2]

        vISweights          = vLL
        # (would need to multiply vLL with prior evaluations and divide by prior evaluations, so latter just cancel)

        vISweightsHelp      = vISweights .- maximum(vISweights)

        vISnormWeights      = exp.(vISweightsHelp) ./ mean(exp.(vISweightsHelp))

        mVarISweights[rr,1] = var(vISnormWeights) #std(filter(isfinite, log.(vISnormWeights))) #


        for pphi = 1:nPhiLast

            vLL                         = vvvaStats[2][pphi][rr][:,2]
            vLLtilde                    = vvvaStats[2][pphi][rr][:,3]

            vISweights                  = vLL - vLLtilde
            # (would need to multiply vLL and vLLtilde both with prior evaluations, so latter just cancel)

            vISweightsHelp              = vISweights .- maximum(vISweights)

            vISnormWeights              = exp.(vISweightsHelp) ./ mean(exp.(vISweightsHelp))

            mVarISweights[rr,pphi+1]    = var(vISnormWeights) #std(filter(isfinite, log.(vISnormWeights))) #

        end

    end

    aVarISweights[:,:,ddgp] = mVarISweights

end



mVarISweightsAllDGPs = reshape(mean(aVarISweights,dims=1), (nPhiLast+1,3))'

mRelVarISweightsAllDGPs = mVarISweightsAllDGPs ./ mVarISweightsAllDGPs[:,1]


plotVarISW       = plot(vXaxis,mVarISweightsAllDGPs[1,:],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="V[w-tilde]",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mVarISweightsAllDGPs[2,:],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mVarISweightsAllDGPs[3,:],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])
xticks!(vXaxis)

savefig(sPlotPath * "plot_VarISW_DGPs.png")



plotRelVarISW       = plot(vXaxis,mRelVarISweightsAllDGPs[1,:],line=(myBlue,0.9,2,vMyShapes[1]),label=vLineLabs[1], xlabel=xLabPhiLast,ylabel="V[w-tilde] relative to LT",title="", xtickfont=font(12),ytickfont=font(12),xguidefontsize=12,yguidefontsize=12,legendfontsize=12,background_color_legend=nothing)

plot!(vXaxis,mRelVarISweightsAllDGPs[2,:],line=(myBlue,0.9,2,vMyShapes[2]),label=vLineLabs[2])
plot!(vXaxis,mRelVarISweightsAllDGPs[3,:],line=(myBlue,0.9,2,vMyShapes[3]),label=vLineLabs[3])
xticks!(vXaxis)

savefig(sPlotPath * "plot_RelVarISW_DGPs.png")
