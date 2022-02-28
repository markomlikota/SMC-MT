# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Jan 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Creates plots of estimation output of one Model specification (DGP-Model-Prior-SMC) over many runs and estimation strategies (LT or MT with different early-stopping parameter values)

# - (marginal) posterior and stage1-density (prior/M0-posterior(w./w.o. early stopping))
# - (marginal) posterior



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec        = "3"     # data from which DGP to take?
SMCspec        = "03"    # which values for SMC tuning parameters to take? (taken to be same for M1 and M0)

vPhiLast       = [0.2, 0.4, 0.6, 0.8, 1.0]

# Remaining options (see main_RBC_OneRun.jl) specified below inside function fRunSMC

nRun            = 1



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


Model1spec      = "L2"
Prior1spec      = "L22"

sPlotPath       = sMyPath * "Output/" * fGetSpecFolder(DGPspec,Model1spec,Prior1spec,SMCspec)

nPhiLast        = length(vPhiLast)

xLabPhiLast     = "\$\\psi_*\$" # "\$\\phi_{N_\\phi}(M_0)\$"

myBlue          = cgrad(:blues)[1.0]
myRed           = cgrad(:reds)[0.8]


include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )




# -- Obtain Time, MDD & #Stages for all estimation specifications --------------


sFileSuffix     = "FinalStats.csv"

vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)



# Time plot:

vTimes          = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"mean")
vTimesStds      = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"std")

plotTimes       = fMyRibbonPlotHere(vPhiLast,vTimes./60,1.64*vTimesStds./60,xLabPhiLast,"minutes")
ylims!((0,ylims(plotTimes)[2]))

savefig(sPlotPath * "plot_meanTimes.png")


vRelTimes       = vTimes/vTimes[1]

plotRelTimes    = fMyRibbonPlotHere(vPhiLast,vRelTimes,1.64*vTimesStds,xLabPhiLast,"time relative to LT")
ylims!((0,ylims(plotRelTimes)[2]))

savefig(sPlotPath * "plot_meanRelTimes.png")



# log MDD plot:

vLogMDDs        = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"mean")
vLogMDDsStds    = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"std")

plotMDDs        = fMyRibbonPlotHere(vPhiLast,vLogMDDs,1.64*vLogMDDsStds,xLabPhiLast)
yLimMean        = mean(ylims(plotMDDs))
ylims!((yLimMean-1,yLimMean+1))

savefig(sPlotPath * "plot_meanLogMDDs.png")


# #  log MDD Std plot:
#
# plotMDDsStds    = fMyRibbonPlotHere(vPhiLast,vLogMDDsStds,zeros(length(vPhiLast)+1),xLabPhiLast)
# ylims!((0,ylims(plotMDDsStds)[2]))
# xticks!([0.0,vPhiLast...])
#
# savefig(sPlotPath * "plot_stdLogMDDs.png")


# #Stages plot:

vNStages        = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,1,"mean",true)

plotNStages     = fMyRibbonPlotHere(vPhiLast,vNStages,zeros(nPhiLast+1),xLabPhiLast)
ylims!((0,ylims(plotNStages)[2]))

savefig(sPlotPath * "plot_meanNStages.png")





# -- Tempering Schedule Plot ---------------------------------------------------


sFileSuffix         = "StageAll_Stats.csv"

vvvaStageAllStats   = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)

plotTSs             = fPlotTSs(vPhiLast,vvvaStageAllStats,"")

savefig(sPlotPath * "plot_medianTSComparison.png")



# -- Posterior Evolution Plots -------------------------------------------------


nRunKDE         = min(nRun,10)


# Object with Posterior Particles:

sFileSuffix     = "StageLast_Particles.csv"
vvvaParticles   = fGetFileAllEstSpecs(DGPspec,SMCspec,nRunKDE,sFileSuffix,vPhiLast)


# Object with Prior Particles:

sFileSuffix     = "Stage1_Particles.csv"
vvvaPrParticles = fGetFileAllEstSpecs(DGPspec,SMCspec,nRunKDE,sFileSuffix,vPhiLast)


# Combine posterior-draws of many runs to form better kernel-density-estimate:

N               = size(vvvaParticles[1][1][1],1) # Number of draws (per run)


aAllParticlesM0 = Array{Float64,3}(undef,N*nRunKDE,nP,nPhiLast+1) #array with M0-posterior-draws in rows (combined from nRunKDE runs), parameters in columns and phiLast (including phiLast=0, so prior) in third dimension

for rr = 1:nRunKDE
    for pphi = 1:nPhiLast
        aAllParticlesM0[(rr-1)*N + 1 : rr*N, :, pphi+1] = vvvaParticles[1][pphi][rr]
    end
    aAllParticlesM0[(rr-1)*N + 1 : rr*N, :, 1] = vvvaPrParticles[1][1][rr]
end

aAllParticlesM1 = Array{Float64,3}(undef,N*nRunKDE,nP,1) #array with all M1-posterior-draws in rows, parameter in columns and phiLast=1 in third dimension
[aAllParticlesM1[(rr-1)*N + 1 : rr*N, :] = vvvaParticles[3][1][rr][:,1:nP] for rr = 1:nRunKDE ]


# Plots:


vBoundTypes = [2,1,2,1,1,1,2,2,1,1,0,1,1,1] # for PrPoPlots; based on parameter domains defined in fIsDrawValid in Priorspec

vDoesParamVaryM1 = [length(unique(aAllParticlesM1[:,pp,1]))!=1 for pp = 1:nP]
vDoesParamVaryM0 = [length(unique(aAllParticlesM0[:,pp,1]))!=1 for pp = 1:nP]
vDoesParamVary = vDoesParamVaryM1 + vDoesParamVaryM0 .== 2

vPIndVary       = collect(1:nP)[vDoesParamVary]


# create vLTXParLabs ...
vLTXParLabs = string.("\$ \\",vParLabs," \$")
vLTXParLabs[2] = "\$ \\sigma_z \$"
vLTXParLabs[4] = "\$ \\sigma_z \$"
vLTXParLabs[5] = "r"
vLTXParLabs[2] = "\$ \\sigma_z \$"
vLTXParLabs[10] = "\$ \\phi_1 \$"
vLTXParLabs[11] = "\$ \\phi_2 \$"

for pInd in vPIndVary

    # Preliminary:

    yLab        = vLTXParLabs[pInd]
    mWeights    = ones(N*nRunKDE,9)
    boundsType  = Int(vBoundTypes[pInd])


    # POSTERIOR EVOLUTION PLOT UNDER M0 (FOR DIFFERENT ϕ_Nϕ_M0) + M1-(FUll-)POSTERIOR + TRUE VALUE (AS 2D PLOT)

    # Posterior evolution under M0 for different phiLast:
    vLineLabs   = string.("\$\\psi_* \$ = ",vPhiLast)

    ppp         = fPlotDistEvol2D(pInd,yLab,aAllParticlesM0[:,:,2:end],mWeights,nPhiLast,boundsType,vLineLabs)
    xLimsToTake = xlims(ppp)
    vLineLabs   = string.("\$\\psi_* \$ = ",[0.0,vPhiLast...])
    ppp         = fPlotDistEvol2D(pInd,yLab,aAllParticlesM0[:,:,1:end],mWeights,nPhiLast+1,boundsType,vLineLabs)
    xlims!(xLimsToTake)

    # Add posterior under M1:
    k           = kde(aAllParticlesM1[:,pInd], weights = Weights(mWeights[:,end]))
    st, en      = fFindPlotBounds(boundsType,k.x)

    if pInd == vPIndVary[end]
        plot!(k.x[st:en], k.density[st:en], line=(myBlue,2.5,:line), label="\$ p(\\theta | Y,M_1) \$", xtickfont=font(14),ytickfont=font(14),ztickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14)
    else
        plot!(k.x[st:en], k.density[st:en], line=(myBlue,2.5,:line), label="\$ p(\\theta | Y,M_1) \$", xtickfont=font(14),ytickfont=font(14),ztickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14,legend=false)
    end

    # # Add red vertical line at true value:
    # vyLimsHere = ylims(ppp)
    # plot!(vec(ones(2,1))*θ0[pInd], [0,vyLimsHere[2]],line=(myRed,2.5,:dot),label="")

    savefig(sPlotPath * string("plot_distEvol2D_p",pInd,".png"))



    ### PLOT OF PRIOR, M0-POSTERIOR (FULL) and M1-(FULL-)POSTERIOR:


    # Plot prior and M0 posterior:
    fPrPoPlot(pInd,aAllParticlesM0,mWeights,boundsType)
    title!(vLTXParLabs[pInd])

    # Add posterior under M1:
    k = kde(aAllParticlesM1[:,pInd], weights = Weights(mWeights[:,end]))
    st, en = fFindPlotBounds(boundsType,k.x)
    plot!(k.x[st:en],k.density[st:en],xlabel=yLab,label="",line=(myBlue,0.9,2,:line),yticks=false)

    # Add red vertical line at true value:
    plot!(vec(ones(length(k.x[st:en]),1))*θ0[pInd],k.density[st:en],label="",line=(myRed,0.9,2,:dot),yticks=false)

    savefig(sPlotPath * string("plot_prpo_p",pInd,".png"))



    ### PLOT OF M0-POSTERIOR (FULL) and M1-(FULL-)POSTERIOR:

    # Plot M1 posterior (still in storage as "k"):
    plot(k.x[st:en],k.density[st:en],xlabel=yLab,label="",line=(myBlue,0.9,2,:line),yticks=false, xtickfont=font(14),ytickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14)

    # Add M0 posterior:
    k = kde(aAllParticlesM0[:,pInd,end], weights = Weights(mWeights[:,end]))
    st, en = fFindPlotBounds(boundsType,k.x)
    plot!(k.x[st:en],k.density[st:en],label="",line=(myBlue,0.9,2,:dash),yticks=false)

    # Add red vertical line at true value:
    plot!(vec(ones(length(k.x[st:en]),1))*θ0[pInd],k.density[st:en],label="",line=(myRed,0.9,2,:dot),yticks=false)

    savefig(sPlotPath * string("plot_popo_p",pInd,".png"))

end




# -- Obtain Importance Sampling Weights for all estimation specifications ------


sFileSuffix     = "Stage1_WLLP.csv"

vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)


mVarISweights   = zeros(nRun,nPhiLast+1)

for rr = 1:nRun


    vLL                 = vvvaStats[3][1][rr][:,2]

    vISweights          = vLL
    # (would need to multiply vLL with prior evaluations and divide by prior evaluations, so latter just cancel)

    vISweightsHelp      = vISweights .- maximum(vISweights)

    vISnormWeights      = exp.(vISweightsHelp) ./ mean(exp.(vISweightsHelp))

    mVarISweights[rr,1] = var(vISnormWeights)


    for pphi = 1:nPhiLast

        vLL                         = vvvaStats[2][pphi][rr][:,2]
        vLLtilde                    = vvvaStats[2][pphi][rr][:,3]

        vISweights                  = vLL - vLLtilde
        # (would need to multiply vLL and vLLtilde both with prior evaluations, so latter just cancel)

        vISweightsHelp              = vISweights .- maximum(vISweights)

        vISnormWeights              = exp.(vISweightsHelp) ./ mean(exp.(vISweightsHelp))

        mVarISweights[rr,pphi+1]    = var(vISnormWeights)

    end

end


vVarISweights = mean(mVarISweights,dims=1)

vRelVarISweights = vVarISweights ./ vVarISweights[1]


plotVarISW      = fMyRibbonPlotHere(vPhiLast,vVarISweights',zeros(nPhiLast+1),xLabPhiLast,"V[w-tilde]")

savefig(sPlotPath * "plot_VarISW.png")


plotRelVarISW       = fMyRibbonPlotHere(vPhiLast,vRelVarISweights',zeros(nPhiLast+1),xLabPhiLast,"V[w-tilde] relative to LT")

savefig(sPlotPath * "plot_RelVarISW.png")
