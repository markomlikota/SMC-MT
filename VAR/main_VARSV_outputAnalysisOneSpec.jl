# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Creates plots of estimation output of one DGP-SMC-specification for all estimation specifications (M1-LT, M1-MT) over many runs.

# Note: as opposed to the file "main_outputAnalysis.jl" in the old folder "SoftwareVAR", here I only consider plots (no tables)



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec         = "3"     # data from which DGP to take?
SMCspec         = "1"     # which values for SMC tuning parameters to take?

vPhiLast        = [0.2, 0.4, 0.6, 0.8, 1.0]

# Remaining options (see main_VARSV_OneRun.jl) specified inside functions to extract files in fFunctionsOutputAnalysis.jl

nRun            = 200



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


using Plots
pyplot()
# using PyPlot
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
sMyPath            = pwd() * "/"


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


Model1spec      = "SV"
Prior1spec      = "SV1"

sPlotPath       = sMyPath * "Output/" * fGetSpecFolder(DGPspec,Model1spec,Prior1spec,SMCspec)

nPhiLast        = length(vPhiLast)

xLabPhiLast     = "\$\\psi_*\$" # "\$\\phi_{N_\\phi}(M_0)\$"

myBlue          = cgrad(:blues)[1.0]
myRed           = cgrad(:reds)[0.8]



# -- Obtain Time, MDD & #Stages for all estimation specifications --------------


sFileSuffix     = "FinalStats.csv"

vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)



### TIME PLOT:


vTimes          = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"mean")
vTimesStds      = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"std")

# plotTimes       = fMyRibbonPlotHere(vPhiLast,vTimes./60,1.64*vTimesStds./60,xLabPhiLast,"minutes")
# ylims!((0,ylims(plotTimes)[2]))

# Rather than assuming symmetric distribution (& in fact normality), compute actual 95% confidence bands:
# (note that the difference is negligible)
vTimesLower     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"lower")
vTimesUpper     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,3,"upper")
vTimesLowerDiff = vTimesLower - vTimes
vTimesUpperDiff = vTimesUpper - vTimes

plotTimes       = fMyRibbonPlotHere(vPhiLast,vTimes./60,(-1*vTimesLowerDiff./60,vTimesUpperDiff./60),xLabPhiLast,"minutes")
ylims!((0,ylims(plotTimes)[2]))

savefig(sPlotPath * "plot_meanTimes.png")


vRelTimes       = vTimes./vTimes[1]
plotRelTimes    = fMyPlot([0,vPhiLast...],vTimes./vTimes[1],"",xLabPhiLast,"time relative to LT","",myBlue)
ylims!((0,1))

savefig(sPlotPath * "plot_meanRelTimes.png")



### LOG MDD PLOT:


vLogMDDs        = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"mean")
vLogMDDsStds    = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"std")

# plotMDDs        = fMyRibbonPlotHere(vPhiLast,vLogMDDs,1.64*vLogMDDsStds,xLabPhiLast)
# yLimMean        = mean(ylims(plotMDDs))
# ylims!((yLimMean-2,yLimMean+2))

vLogMDDsLower     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"lower")
vLogMDDsUpper     = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,4,"upper")
vLogMDDsLowerDiff = vLogMDDsLower - vLogMDDs
vLogMDDsUpperDiff = vLogMDDsUpper - vLogMDDs

plotMDDs        = fMyRibbonPlotHere(vPhiLast,vLogMDDs,(-1*vLogMDDsLowerDiff,vLogMDDsUpperDiff),xLabPhiLast)
yLimMean        = mean(ylims(plotMDDs))
ylims!((yLimMean-2,yLimMean+2))

savefig(sPlotPath * "plot_meanLogMDDs.png")



###  STD LOG MDD PLOT:


plotMDDsStds    = fMyRibbonPlotHere(vPhiLast,vLogMDDsStds,zeros(length(vPhiLast)+1),xLabPhiLast)
ylims!((0,ylims(plotMDDsStds)[2]))
xticks!([0.0,vPhiLast...])

savefig(sPlotPath * "plot_stdLogMDDs.png")



### #STAGES PLOT:


vNStages        = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,1,"mean",true)
vNStagesLower   = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,1,"lower",true)
vNStagesUpper   = fGetStatIn_vvvaObject(vvvaStats,vPhiLast,1,"upper",true)

plotNStages     = fMyRibbonPlotHere(vPhiLast,vNStages,(vNStages-vNStagesLower,vNStagesUpper-vNStages),xLabPhiLast)
ylims!((0,ylims(plotNStages)[2]))

savefig(sPlotPath * "plot_meanNStages.png")



# -- Tempering Schedule Plot ---------------------------------------------------


sFileSuffix         = "StageAll_Stats.csv"

vvvaStageAllStats   = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)

plotTSs             = fPlotTSs(vPhiLast,vvvaStageAllStats,"")

savefig(sPlotPath * "plot_medianTSComparison.png")



# -- Parameter Posterior Moment Plots ------------------------------------------


sFileSuffix         = "StageLast_ParticleStats.csv"

vvvaStats           = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)


# Get true parameters:

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )
include(sMyPath * "SpecFiles/" * string("script_ModelspecSV.jl") )
θ0          = fVARObjectsToϑ_SV(Φ0, Σ0, vρ0, vξ0)


# Colors and line shapes:

mColsShapes             = Matrix{Any}(undef,13,2)

mColsShapes[mParaPointers_SV[1,1]:(Int(mParaPointers_SV[1,2]/2)),1]   .= cgrad(:grays)[0.9]
mColsShapes[(Int(mParaPointers_SV[1,2]/2)+1):mParaPointers_SV[1,2],1] .= cgrad(:grays)[0.9]
mColsShapes[mParaPointers_SV[2,1]:mParaPointers_SV[2,2],1]            .= cgrad(:grays)[0.7]

mColsShapes[mParaPointers_SV[3,1]:mParaPointers_SV[3,2],1]            .= cgrad(:blues)[0.5]
mColsShapes[mParaPointers_SV[4,1]:mParaPointers_SV[4,2],1]            .= cgrad(:blues)[1.0]

myShapes                = [:solid, :dash, :dot]

mColsShapes[1:9,2]      = repeat([:solid],9,1)#repeat(myShapes,3,1)
mColsShapes[10:13,2]    = repeat(myShapes[1:2],2,1)

vThickness              = [repeat([2],9,1);repeat([3],4,1)]


# Latex Labels :

mLTXLabs_Φ1    = string.(repeat(["Phi_{1,"],n,n) .* repeat(string.(collect(1:n)),1,n) .* repeat(string.(collect(1:n)'),n,1) .* "}")
mLTXLabs_Φc    = string.(repeat(["Phi_{c,"],n,1) .* string.(collect(1:n)) .* "}")
mLTXLabs_Φ     = permutedims([mLTXLabs_Φ1 mLTXLabs_Φc])

mLTXLabs_Σ    = string.(repeat(["Sigma_{"],n,n) .* repeat(string.(collect(1:n)),1,n) .* repeat(string.(collect(1:n)'),n,1) .* "}")

vLTXLabs_vρ    = string.(repeat(["rho_"],n,1) .* string.(collect(1:n)))

vLTXLabs_vξ    = string.(repeat(["xi_"],n,1) .* string.(collect(1:n)))

vLTXParLabs    = [vec(mLTXLabs_Φ)..., mLTXLabs_Σ[tril(trues(size(mLTXLabs_Σ)))]..., vLTXLabs_vρ..., vLTXLabs_vξ...]

vLTXParLabs = string.("\$\\",vLTXParLabs,"\$")


# Plots:

function fPlotParticlePostStat(vpp,ss,vvvaStats,sCaption)

    pp=vpp[1]
    vPostStat = [mean(getindex.(vvvaStats[3][1],pp,ss)),[mean(getindex.(vvvaStats[2][pphi],pp,ss)) for pphi = 1:nPhiLast]...]
    ssPlot = plot([0,vPhiLast...],vPostStat,line=(mColsShapes[pp,1],0.9,vThickness[pp],mColsShapes[pp,2]),label=vLTXParLabs[pp], legend = :outerright,xlabel=xLabPhiLast,title=sCaption,xtickfont=font(14),ytickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=11,background_color_legend = nothing)
    xticks!([0.0,vPhiLast...])
    for pp = vpp[2:end]
        vPostStat = [mean(getindex.(vvvaStats[3][1],pp,ss)),[mean(getindex.(vvvaStats[2][pphi],pp,ss)) for pphi = 1:nPhiLast]...]
        plot!([0,vPhiLast...],vPostStat,line=(mColsShapes[pp,1],0.9,vThickness[pp],mColsShapes[pp,2]),label=vLTXParLabs[pp])
    end

    return ssPlot

end


plotPostMeans = fPlotParticlePostStat(1:13,1,vvvaStats,"Posterior Means")
savefig(sPlotPath * "plot_post_mean.png")

plotPostVars = fPlotParticlePostStat(1:13,2,vvvaStats,"Posterior Variance")
savefig(sPlotPath * "plot_post_var.png")

plotPostLower = fPlotParticlePostStat(1:13,3,vvvaStats,"Posterior 5th percentile")
savefig(sPlotPath * "plot_post_lower.png")

plotPostUpper = fPlotParticlePostStat(1:13,4,vvvaStats,"Posterior 95th percentile")
savefig(sPlotPath * "plot_post_upper.png")



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


aAllParticlesM0 = Array{Float64,3}(undef,N*nRunKDE,9,nPhiLast+1) #array with M0-posterior-draws in rows (combined from nRunKDE runs), parameters in columns and phiLast (including phiLast=0, so prior) in third dimension

for rr = 1:nRunKDE
    for pphi = 1:nPhiLast
        aAllParticlesM0[(rr-1)*N + 1 : rr*N, :, pphi+1] = vvvaParticles[1][pphi][rr]
    end
    aAllParticlesM0[(rr-1)*N + 1 : rr*N, :, 1] = vvvaPrParticles[1][1][rr]
end

aAllParticlesM1 = Array{Float64,3}(undef,N*nRunKDE,9,1) #array with all M1-posterior-draws in rows, parameter in columns and phiLast=1 in third dimension
[aAllParticlesM1[(rr-1)*N + 1 : rr*N, :] = vvvaParticles[3][1][rr][:,1:9] for rr = 1:nRunKDE ]


# Plots:

vBoundTypes     = zeros(9)
vBoundTypes[[7,9]] .=1
global vParLabs = vLTXParLabs

for pInd = 1:9

    # Preliminary:

    yLab        = vLTXParLabs[pInd]
    mWeights    = ones(N*nRunKDE,9)
    nLines      = 5#6
    boundsType  = Int(vBoundTypes[pInd])



    ### POSTERIOR EVOLUTION PLOT UNDER M0 (FOR DIFFERENT ϕ_Nϕ_M0) + M1-(FUll-)POSTERIOR + TRUE VALUE:


    # Posterior evolution under M0 for different phiLast:
    vLineLabs   = string.("\$\\psi_* \$ = ",vPhiLast)

    ppp         = fPlotDistEvol2D(pInd,yLab,aAllParticlesM0[:,:,2:6],mWeights,nLines,boundsType,vLineLabs)
    xLimsToTake = xlims(ppp)
    vLineLabs   = string.("\$\\psi_* \$ = ",[0.0,vPhiLast...])
    ppp         = fPlotDistEvol2D(pInd,yLab,aAllParticlesM0[:,:,1:6],mWeights,6,boundsType,vLineLabs)
    xlims!(xLimsToTake)

    # Add posterior under M1:
    k           = kde(aAllParticlesM1[:,pInd], weights = Weights(mWeights[:,end]))
    st, en      = fFindPlotBounds(boundsType,k.x)

    if pInd == 9
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
    plot!(k.x[st:en],k.density[st:en],label="",line=(myBlue,0.9,2,:line),yticks=false)

    # Add red vertical line at true value:
    plot!(vec(ones(length(k.x[st:en]),1))*θ0[pInd],k.density[st:en],label="",line=(myRed,0.9,2,:dot),yticks=false)

    savefig(sPlotPath * string("plot_prpo_p",pInd,".png"))



    ### PLOT OF M0-POSTERIOR (FULL) and M1-(FULL-)POSTERIOR:

    # Plot M1 posterior (still in storage as "k"):
    plot(k.x[st:en],k.density[st:en],label="",line=(myBlue,0.9,2,:line),yticks=false, xtickfont=font(14),ytickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14)

    # Add M0 posterior:
    k = kde(aAllParticlesM0[:,pInd,end], weights = Weights(mWeights[:,end]))
    st, en = fFindPlotBounds(boundsType,k.x)
    plot!(k.x[st:en],k.density[st:en],label="",line=(myBlue,0.9,2,:dash),yticks=false)

    # Add red vertical line at true value:
    plot!(vec(ones(length(k.x[st:en]),1))*θ0[pInd],k.density[st:en],label="",line=(myRed,0.9,2,:dot),yticks=false)

    savefig(sPlotPath * string("plot_popo_p",pInd,".png"))

end



# -- Importance Sampling Weights -----------------------------------------------


sFileSuffix     = "Stage1_WLLP.csv"

vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)


# Obtain variance for each run:

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


# For each specification compute average variance across runs, & 5th and 95h percentiles:

vVarISweights = mean(mVarISweights,dims=1)

mSortedVarISweights = sort(mVarISweights,dims=1)

vVarISweights_lower = mSortedVarISweights[Int(nRun*0.05),:]
vVarISweights_upper = mSortedVarISweights[Int(nRun*0.95),:]

vVarISweights_lowerDiff = vVarISweights_lower - vVarISweights'
vVarISweights_upperDiff = vVarISweights_upper - vVarISweights'


plotVarISW      = fMyRibbonPlotHere(vPhiLast,vVarISweights',(-1*vVarISweights_lowerDiff,vVarISweights_upperDiff),xLabPhiLast,"V[w-tilde]")
ylims!((0,ylims(plotVarISW)[2]))
#plotVarISW      = fMyRibbonPlotHere(vPhiLast,vVarISweights',zeros(nPhiLast+1),xLabPhiLast,"V[w-tilde]")

savefig(sPlotPath * "plot_VarISW.png")


# Relative to LT:

vRelVarISweights = vVarISweights ./ vVarISweights[1]

plotRelVarISW       = fMyRibbonPlotHere(vPhiLast,vRelVarISweights',zeros(nPhiLast+1),xLabPhiLast,"V[w-tilde] relative to LT")
ylims!((0,ylims(plotRelVarISW)[2]))
savefig(sPlotPath * "plot_RelVarISW.png")



# -- Analyze changes in benefit under changed runtime --------------------------



### #Stages for M0 and M1 as functions of phi:

sFileSuffix     = "FinalStats.csv"

vvvaStats       = fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)


mNtilde         = zeros(nPhiLast+1,2)

mNtilde[:,2]      = mean.( [getindex.(vvvaStats[3][1],1),[getindex.(vvvaStats[2][pphi],1) for pphi = 1:nPhiLast]...] ) .+ 1
mNtilde[2:end,1]  = mean.( [getindex.(vvvaStats[1][pphi],1) for pphi = 1:nPhiLast] ) .+ 1


# Plot:

vXaxis = [0,vPhiLast...]


fMyPlot(vXaxis,mNtilde[:,1],"",xLabPhiLast,"","\$ \\tilde{N}_0(\\psi_*) \$",myRed)
plot!(vXaxis,mNtilde[:,2],line=(myBlue,0.9,2,:line),label="\$ \\tilde{N}_1(\\psi_*) \$")
xticks!(vXaxis)

savefig(sPlotPath * "plot_NStagesM0M1.png")



### Analytically approximated runtimes:


# Obtain LL-times (see main_VARSV_LLevaltimes.jl):

tM0     = 0.005460195
tM1     = 0.049899417
t       = 9.13876


# Function to compute approximated runtime reduction:

function fTimesM1(tM1,tM0)

    vTimesM1 = 500 * [mNtilde[1,2]*tM1, [mNtilde[pphi,1] * tM0 + mNtilde[pphi,2] * (tM1 + tM0) for pphi = 2:nPhiLast+1]...]

    return vTimesM1./vTimesM1[1]

end


# Plot:

vFactor = [2,4,1e100]

M = length(vFactor)+1
mstep = 1
vBlues = [cgrad(:blues)[mm/(M*mstep)] for mm = 1:mstep:M*mstep]

sRelTs = "\$ \\tau_1/\\tau_0 \$ = "
fRepresentNumber(tt) = Int(floor(tt))

relTaus1 = tM1/tM0
yLab = "" # "\$ \\mathcal{R}(\\psi_*,\\tau_0/\\tau_1) \$"
fMyPlot(vXaxis,fTimesM1(tM1,tM0),"",xLabPhiLast,yLab,string(sRelTs,fRepresentNumber(relTaus1)),vBlues[1])

plot!(vXaxis,fTimesM1(tM1*vFactor[1],tM0),line=(vBlues[1+1],0.9,2,:line),label=string(sRelTs,fRepresentNumber(relTaus1*vFactor[1])))
plot!(vXaxis,fTimesM1(tM1*vFactor[2],tM0),line=(vBlues[1+2],0.9,2,:line),label=string(sRelTs,fRepresentNumber(relTaus1*vFactor[2])))
plot!(vXaxis,fTimesM1(tM1*vFactor[3],tM0),line=(vBlues[1+3],0.9,2,:line),label="\$ \\tau_1/\\tau_0 \\rightarrow \\infty \$")
xticks!(vXaxis)
ylims!((0,1))


savefig(sPlotPath * "plot_relTimes_Approx.png")



# Notes on comparing approximated runtimes (not reductions) to actual:

# Analytically approximated runtime underestimates actual runtime for lower phis...
# ... (and hence underestimates runtime reduction due to MT)

# Adjusting for the fact that some draws are invalid in first stage doesn't help.
# So underestimation is likely due to the LL-eval being slower for particles that have a very low LL, and there are more such particles in initial stages...
