# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Tests SMC function on simple example: p(Y|θ) = N(θ,1), p(θ) = N(0,1)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


nProcs = 2

using Distributed
addprocs(nProcs-length(procs()))


using Plots
pyplot()
using StatsBase
using DataFrames
using CSV

@everywhere using ParallelDataTransfer
@everywhere using SharedArrays
@everywhere using Distributions
@everywhere using LinearAlgebra



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


@everywhere cd()
@everywhere sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareVAR/")


include(sMyPath * "/Functions/fSetupSMC.jl")
@everywhere include(sMyPath * "/Functions/resampleNYFED.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


# Specify path output:

sSpecFolder     = "SMCtest"
try
    mkdir(sMyPath * "Output/" * sSpecFolder)
catch
end

sSpecPath = sMyPath * "Output/" * sSpecFolder * "/"


# Specify functions for SMC:

y       = 0.4

@eval @everywhere y = $y

@everywhere begin
    fLLtilde(ϑ)         = 0
    fLL(ϑ)              = log(pdf(Normal(ϑ[1],1),y))
    fPriorDraw()        = rand(Normal())
    fIsDrawValid(ϑ)     = true
    fPriorLogEval(ϑ)    = log(pdf(Normal(),ϑ[1]))
end


# For computation of true MDD (for plots):

fTrueLogMDD(y)  = -log(2π)/2 - log(2)/2 - 1/4 * y^2


# SMC settings:

N               = 500
useAdaptive     = 0
λ               = 2
Nϕ              = 5
ϕ_Nϕ            = 1.00 # value for phi in last stage (fixed TS can target it exactly; adaptive TS terminates when phi becomes larger or equal to ϕ_Nϕ)
α               = 0.95
Nmh             = 1
nB              = 1 # number of blocks for RWMH in mutation steps, random blocking is used
nP              = 1
c0              = 0.5
accStar         = 0.25
Nbar            = N/2
showMessages    = false

tSettings       = (N,useAdaptive,λ,Nϕ,ϕ_Nϕ,α,Nmh,nB,nP,c0,accStar,Nbar,showMessages)
@everywhere @eval tSettings = $tSettings
vParLabs        = ["theta"]


sOutputPath         = sSpecPath
sOutputFilePrefix   = ""
setSeed             = 0


# Run SMC many times:

M           = 1000

vLogMDD     = zeros(M,1)
vPostMean   = zeros(M,1)
vPostVar    = zeros(M,1)

for mm = 1:M

    println(string("SMC Test Run ",Int(mm)))

    aParticles, mWeights, aLogLiks, mPriors, vRejRate, vSstep, vESS, vc, vϕ, logMDD, timeTotal = fSMC(tSettings,sOutputPath,sOutputFilePrefix,setSeed)

    vLogMDD[mm]     = sum(logMDD)
    vPostMean[mm]   = mean(aParticles[:,1,end])
    vPostVar[mm]    = var(aParticles[:,1,end])

end


# Plots:

histogram(vLogMDD,label="",title="SMC-implied logMDD vs true value")
vline!([fTrueLogMDD(y)],line =(:red,0.9,3),label="")
savefig(sSpecPath * "plot_mdd.png")

histogram(vPostMean,label="",title="SMC-implied posterior mean vs true value")
vline!([y/2],line =(:red,0.9,3),label="")
savefig(sSpecPath * "plot_postmean.png")

histogram(vPostVar,label="",title="SMC-implied posterior variance vs true value")
vline!([1/2],line =(:red,0.9,3),label="")
savefig(sSpecPath * "plot_postvar.png")
