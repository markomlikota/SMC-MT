# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Main script to simulate data for RBC based on model "L1"
# (first order linearization of RBC with cap. adj. costs),
# given DGP specified below.


# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !
# i.e. change definitions of fMyCSVREAD, fMyCSVWRITE and fMyCSVAPPEND in "fFolderFileManagement.jl"



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec     = "1"         # which DGP to simulate?

writeOutput = 1           # 1: yes, 0: only simulate, do not store neither data as csv nor plots
useModelspecInFileName = 1 # KEEP HERE AT DEFAULT: 1 !!!!
# 1: yes, 0: no (only data without Modelspec in filename can be used for estimation)

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



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


cd()
sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/RBC/")
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
    σy_0  = 0
    σi_0  = 0
    σl_0  = 0
    θ0          = [ρz_0, σz_0, ρb_0, σb_0, r_0, τ_0, α_0, δ_0, ν_0, ϕ1_0, ϕ2_0, σy_0, σi_0, σl_0]
end


# Load Model specification file (based on which simulation is performed):

Modelspec       = "L1"
include(sMyPath * "SpecFiles/" * string("script_Modelspec",Modelspec,".jl") )


# Simulate data:

Ψ0, Ψ1, Ψu, Φ1, Φε, vYlabs, vXlabs  = fGetStateSpaceRep(θ0)

Random.seed!(354)

mY, mX, mε, mu      = fSimulateLinearSS(Ψ0, Ψ1, Ψu, Φ1, Φε, T)


# Store data (measurement variables y):

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

p1          = fMyPlot(1:T,mY[:,1],vYlabs[1],"","")
p2          = fMyPlot(1:T,mY[:,2],vYlabs[2],"","")
p3          = fMyPlot(1:T,mY[:,3],vYlabs[3],"t","")
pp          = plot(p1,p2,p3,layout=(3,1),legend=false)
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


# Plot all (econometric) state variables individually:

vPlots = Array{Any}(undef,1,size(mX,2))

for xx in 1:size(mX,2)

    xx == 6 || xx == 7 ? syLabHere = "" : syLabHere = "% dev from SS"
    xx == 6 || xx == 7 ? mData = mX[:,xx] : mData = 100*mX[:,xx]

    vPlots[xx]                 = fMyPlot(0:T,mData,vXlabs[xx],"t",syLabHere)

    if writeOutput == 1

        if useModelspecInFileName == 0
            sFileForOutput         = "Data_" * string("dgp",DGPspec) * string("_var_",vXlabs[xx]) * ".png"
        else
            sFileForOutput         = "Data_" * string("dgp",DGPspec) * string("_m",Modelspec)* string("_var_",vXlabs[xx]) * ".png"
        end
        sPathForOutput         = sFolderForOutput * sFileForOutput

        savefig(sPathForOutput)

    end

end


# Compute sample variance of y^o, i^o, l^o to determine appropriate measurement error variances:

if noMEs # (only done in case DGP without measurement errors is used in this run)

    vSampleVars     = var(mY,dims=1)

    vMEVarToUse     = (fracVarMEs/(1-fracVarMEs)) * vSampleVars

    vMEStDsToUse    = sqrt.(vMEVarToUse)

    println(string("Use σy_0 = ",round(vMEStDsToUse[1];digits=16)))
    println(string("Use σi_0 = ",round(vMEStDsToUse[2];digits=16)))
    println(string("Use σl_0 = ",round(vMEStDsToUse[3];digits=16)))

else

    vSampleVars     = var(mY,dims=1)

    println(string("σy^2 is ",round( 100* σy_0^2 / vSampleVars[1];digits=4),"% of total sample variance of ln y^o"))
    println(string("σi^2 is ",round( 100* σi_0^2 / vSampleVars[2];digits=4),"% of total sample variance of ln i^o"))
    println(string("σl^2 is ",round( 100* σl_0^2 / vSampleVars[3];digits=4),"% of total sample variance of ln l^o"))

end
