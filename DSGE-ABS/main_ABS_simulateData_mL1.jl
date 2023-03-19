# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# May 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Main script to simulate data for ABS based on model "L1" given DGP specified below.


# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !
# i.e. change definitions of fMyCSVREAD, fMyCSVWRITE and fMyCSVAPPEND in "fFolderFileManagement.jl"



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec     = "0"         # which DGP to simulate?

writeOutput = 1           # 1: yes, 0: only simulate, do not store neither data as csv nor plots
useModelspecInFileName = 1 # KEEP HERE AT DEFAULT: 1 !!!!
# 1: yes, 0: no (only data without Modelspec in filename can be used for estimation)

noMEs       = false        # set measurement errors to zero? -> e.g. for calibrating them
fracVarMEs  = 0.10        # (in case the script is used to calibrate StDs of MEs: which fraction of sample variance of observables should be due to MEs?)


randomSeedNumber    = 349
randomSeedNumber2   = 24


# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


using Plots
pyplot()
using DataFrames
using CSV

using Random
using LinearAlgebra
using Statistics



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


cd()
sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/DSGE-ABS/")
# sMyPath            = pwd() * "/"


include(sMyPath * "/Functions/fSimulateLinearSS.jl")
include(sMyPath * "/Functions/fMyPlot.jl")
include(sMyPath * "/Functions/fFolderFileManagement.jl")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


sFolderForOutput         = sMyPath * "Data/"


# Load DGP specification file:

include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )

if noMEs
    σ_dy_0       = 0
    σ_dwnom_0    = 0
    σ_r_0        = 0
    σ_infl_0     = 0
    θ0          = [rA_0, πA_0, γA_0, τ_0, ν_0, κ_0, φw_0, ψw_0, ψp_0, ψ1_0, ψ2_0, ρr_0, ρg_0, ρa_0, ρp_0, σr_0, σg_0, σa_0, σp_0, gStar_0, λpSS_0, λwSS_0, χh_0, σ_dy_0, σ_dwnom_0, σ_r_0, σ_infl_0]
end


# Load Model specification file (based on which simulation is performed):

Modelspec       = "L1"
include(sMyPath * "SpecFiles/" * string("script_Modelspec",Modelspec,".jl") )


# Simulate data:

mY = fSimulate_mL1(θ0,T,randomSeedNumber,randomSeedNumber2)
# Ψ0, Ψ1, Φ1, Φε, Σ_u = fSolveModelL1(θ0)
#
# Ψu = diagm([σ_dy_0, σ_dwnom_0, σ_r_0, σ_infl_0])
#
# Random.seed!(354)
#
# mY, mX, mε, mu      = fSimulateLinearSS(Ψ0, Ψ1, Ψu, Φ1, Φε, T)


# Store data (measurement variables y):

vYlabs      = ["YGR^o","WGR^o","INT^o","INFL^o"]

if writeOutput == 1

    if useModelspecInFileName == 0
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * ".csv"
    else
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * string("_m",Modelspec) * ".csv"
    end
    sPathForOutput         = sFolderForOutput * sFileForOutput

    # CSV.write(sPathForOutput,  DataFrame(permutedims(vYlabs)), header=false)
    # CSV.write(sPathForOutput,  DataFrame(mY), append=true)
    fMyCSVWRITE(sPathForOutput,mY,vYlabs)

end


# Plot data (measurement variables y):

vYlatexLabs = string.("\$ \\; ",["YGR", "WGR", "INT", "INFL"],"^o \$")

myBlue      = cgrad(:blues)[1.0]

p1          = fMyPlot(1:T,mY[:,1],vYlatexLabs[1],"","","",myBlue)
# yticks!([-0.025,0.025,0.075])
p2          = fMyPlot(1:T,mY[:,2],vYlatexLabs[2],"","","",myBlue)
# yticks!([-1.16,-1.13,-1.10])
p3          = fMyPlot(1:T,mY[:,3],vYlatexLabs[3],"t","","",myBlue)
# yticks!([-0.04,0.00,0.04])
p4          = fMyPlot(1:T,mY[:,4],vYlatexLabs[4],"t","","",myBlue)
pp          = plot(p1,p2,p3,p4,layout=(2,2),legend=false)
display(pp)

if writeOutput == 1

    if useModelspecInFileName == 0
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * ".png"
    else
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * string("_m",Modelspec) * ".png"
    end
    sPathForOutput         = sFolderForOutput * sFileForOutput

    savefig(sPathForOutput)

end



# Compute sample variance of observables to determine appropriate measurement error variances:

if noMEs # (only done in case DGP without measurement errors is used in this run)

    vSampleVars     = var(mY,dims=1)

    vMEVarToUse     = (fracVarMEs/(1-fracVarMEs)) * vSampleVars

    vMEStDsToUse    = sqrt.(vMEVarToUse)

    println("Use:")
    println(string("σ_dy_0 = ",round(vMEStDsToUse[1];digits=16)))
    println(string("σ_dwnom_0 = ",round(vMEStDsToUse[2];digits=16)))
    println(string("σ_r_0 = ",round(vMEStDsToUse[3];digits=16)))
    println(string("σ_infl_0 = ",round(vMEStDsToUse[4];digits=16)))

else

    vSampleVars     = var(mY,dims=1)

    println(string("σ_dy^2 is ",round( 100* σ_dy_0^2 / vSampleVars[1];digits=4),"% of total sample variance of YGR^o"))
    println(string("σ_dwnom^2 is ",round( 100* σ_dwnom_0^2 / vSampleVars[2];digits=4),"% of total sample variance of WGR^o"))
    println(string("σ_r^2 is ",round( 100* σ_r_0^2 / vSampleVars[3];digits=4),"% of total sample variance of INT^o"))
    println(string("σ_infl^2 is ",round( 100* σ_infl_0^2 / vSampleVars[4];digits=4),"% of total sample variance of INFL^o"))

end
