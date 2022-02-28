# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Creates plots of estimation output of several specifications (DGP-Model-Prior-SMC) for LT and MT.



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


# HARDCODED (FOR NOW)
# vDGPspecs         = ["1","2","3"]       # data from which DGP to take?



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
sMyPath            = string(pwd(),"/Dropbox/ModelTempering/Software/RBC/")
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


sFolderForOutput       = sMyPath * "Data/"

myBlue          = cgrad(:blues)[1.0]
myRed           = cgrad(:reds)[0.8]
myGreen         = cgrad(:greens)[0.8]


### Compare data from multiple DGPs:

DGPspec         = 1
include(sMyPath * "SpecFiles/" * string("script_DGPspec",DGPspec,".jl") )


aData   = Array{Float64}(undef,T,3,3)

[ aData[:,:,dd] = fGetFile(sFolderForOutput, "Data_" * string("dgp",dd) * ".csv") for dd = 1:3 ]


f = open(sFolderForOutput * "Data_dgp1.csv")
vYlabs = readdlm(f, ',',skipstart=0)[1,:]


p1          = fMyPlot(1:T,aData[:,1,1],vYlabs[1],"","","DGP 1",myBlue)
plot!(1:T,aData[:,1,2],label="DGP 2",line=(myRed,0.9,2,:line))
plot!(1:T,aData[:,1,3],label="DGP 3",line=(myGreen,0.9,2,:dash))

p2          = fMyPlot(1:T,aData[:,2,1],vYlabs[2],"","","DGP 1",myBlue)
plot!(1:T,aData[:,2,2],label="DGP 2",line=(myRed,0.9,2,:line))
plot!(1:T,aData[:,2,3],label="DGP 3",line=(myGreen,0.9,2,:dash))

p3          = fMyPlot(1:T,aData[:,3,1],vYlabs[3],"","","DGP 1",myBlue)
plot!(1:T,aData[:,3,2],label="DGP 2",line=(myRed,0.9,2,:line))
plot!(1:T,aData[:,3,3],label="DGP 3",line=(myGreen,0.9,2,:dash))

# ,extra_kwargs=:subplot,legend_position = :outerbottom,legend_column=-1

pp          = plot(p1,p2,p3,layout=(3,1),legend=:none)
display(pp)

sFileForOutput         = "Data_" * "allDGPs" * ".png"
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)
