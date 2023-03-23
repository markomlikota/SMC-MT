# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Jan 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Main script to simulate data for RBC based on model "L2" given DGP specified below.


# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !
# i.e. change definitions of fMyCSVREAD, fMyCSVWRITE and fMyCSVAPPEND in "fFolderFileManagement.jl"



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec     = "3"         # which DGP to simulate?

writeOutput = 1           # 1: yes, 0: only simulate, do not store neither data as csv nor plots
useModelspecInFileName = 2 # BEST TO KEEP AT DEFAULT; 2
# 1: yes, 0: no (only data without Modelspec in filename can be used for estimation), 2: do both

noMEs       = false        # set measurement errors to zero? -> e.g. for calibrating them
fracVarMEs  = 0.05        # (in case the script is used to calibrate StDs of MEs: which fraction of sample variance of observables should be due to MEs?)


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
using Distributions
using Interpolations
using Roots
using SolveDSGE



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


cd()
sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/RBC/")
# sMyPath            = pwd() * "/"


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
    σy_0  = 0
    σi_0  = 0
    σl_0  = 0
    θ0          = [ρz_0, σz_0, ρb_0, σb_0, r_0, τ_0, α_0, δ_0, ν_0, ϕ1_0, ϕ2_0, σy_0, σi_0, σl_0]
end


# Load Model specification file (based on which simulation is performed):

Modelspec       = "L2"
include(sMyPath * "SpecFiles/" * string("script_Modelspec",Modelspec,".jl") )


# Simulate data:

fΦ, fΨ, fDrawInitialStatesBSPF = fSolveModelL2(θ0)

nY      = 3
nS      = 3
nε      = 2

Random.seed!(354)

mY, mX, mε, mu = fSimulate_mL2(fΨ,fΦ,θ0,nY,nS,nε, T)


# Store data (measurement variables y):

vYlabs      = ["log y^o","log i^o","log l^o"]

if writeOutput == 1
    if useModelspecInFileName == 0 || useModelspecInFileName == 2
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * ".csv"
        sPathForOutput         = sFolderForOutput * sFileForOutput

        fMyCSVWRITE(sPathForOutput,mY,vYlabs)
    end
    if useModelspecInFileName == 1 || useModelspecInFileName == 2
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * string("_m",Modelspec) * ".csv"
        sPathForOutput         = sFolderForOutput * sFileForOutput

        fMyCSVWRITE(sPathForOutput,mY,vYlabs)
    end

end


# Plot data (measurement variables y):

vYlatexLabs = string.("\$ ln \\; ",["Y", "I", "L"],"^o \$")

myBlue      = cgrad(:blues)[1.0]

p1          = fMyPlot(1:T,mY[:,1],vYlatexLabs[1],"","","",myBlue)
# yticks!([-0.025,0.025,0.075])
p2          = fMyPlot(1:T,mY[:,2],vYlatexLabs[2],"","","",myBlue)
# yticks!([-1.16,-1.13,-1.10])
p3          = fMyPlot(1:T,mY[:,3],vYlatexLabs[3],"t","","",myBlue)
# yticks!([-0.04,0.00,0.04])
pp          = plot(p1,p2,p3,layout=(3,1),legend=false)
display(pp)

cSS, kSS, lSS, ySS, iSS, Zstar, Bstar = fGetSteadyState(θ0)

if writeOutput == 1

    if useModelspecInFileName == 0 || useModelspecInFileName == 2
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * ".png"
        sPathForOutput         = sFolderForOutput * sFileForOutput
        savefig(sPathForOutput)
    end
    if useModelspecInFileName == 1 || useModelspecInFileName == 2
        sFileForOutput         = "Data_" * string("dgp",DGPspec) * string("_m",Modelspec) * ".png"
        sPathForOutput         = sFolderForOutput * sFileForOutput
        savefig(sPathForOutput)
    end

end



# Compute sample variance of y^o, i^o, l^o to determine appropriate measurement error variances:

if noMEs # (only done in case DGP without measurement errors is used in this run)

    vSampleVars     = var(mY,dims=1)

    vMEVarToUse     = (fracVarMEs/(1-fracVarMEs)) * vSampleVars

    vMEStDsToUse    = sqrt.(vMEVarToUse)

    println("Use:")
    println(string("σy_0 = ",round(vMEStDsToUse[1];digits=16)))
    println(string("σi_0 = ",round(vMEStDsToUse[2];digits=16)))
    println(string("σl_0 = ",round(vMEStDsToUse[3];digits=16)))

else

    vSampleVars     = var(mY,dims=1)

    println(string("σy^2 is ",round( 100* σy_0^2 / vSampleVars[1];digits=4),"% of total sample variance of ln y^o"))
    println(string("σi^2 is ",round( 100* σi_0^2 / vSampleVars[2];digits=4),"% of total sample variance of ln i^o"))
    println(string("σl^2 is ",round( 100* σl_0^2 / vSampleVars[3];digits=4),"% of total sample variance of ln l^o"))

end
