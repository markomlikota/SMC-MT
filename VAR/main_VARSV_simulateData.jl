# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Main script to simulate data for VAR-SV based on assumed DGP.



# -------------------------------------------------------------------------------

# OPTIONS TO SPECIFY:
# -------------------------------------------------------------------------------


DGPspec   = "3"         # which DGP to simulate?




# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


using Plots
pyplot()
using LaTeXStrings
using DataFrames
using CSV

using Random
using Distributions
using LinearAlgebra



# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


# cd()
# sMyPath            = string(pwd(),"/Dropbox/FileExchange_FS_MM/SMC-MT/SoftwareVAR/")
sMyPath            = pwd() * "/"


include(sMyPath * "/Functions/fSimulateVARSVp1.jl")
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


# Simulate data:

Random.seed!(3254)
mData, vmΣ, md  = fSimulateVARSVp1(Φ0,Φε0,vd0,vρ0,vξ0,T)


sFileForOutput         = "Data_" * string("dgp",DGPspec) * ".csv"
sPathForOutput         = sFolderForOutput * sFileForOutput

fMyCSVWRITE(sPathForOutput,mData,["y1","y2"])


sFileForOutput         = "SV_" * string("dgp",DGPspec) * ".csv"
sPathForOutput         = sFolderForOutput * sFileForOutput

fMyCSVWRITE(sPathForOutput,md',["d1","d2"])


# Plot y:

p1          = fMyPlot(1:T,mData[:,1],"\$y_{1t}\$","")
p2          = fMyPlot(1:T,mData[:,2],"\$y_{2t}\$","t")
pp          = plot(p1,p2,layout=(2,1),legend=false)
display(pp)


sFileForOutput         = "Data_" * string("dgp",DGPspec) * ".png"
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)


# Plot d:
p1          = fMyPlot(1:T,md'[:,1],"\$d_{1t}\$","")
p2          = fMyPlot(1:T,md'[:,2],"\$d_{2t}\$","t")
pp          = plot(p1,p2,layout=(2,1),legend=false)
display(pp)

sFileForOutput         = "SV_" * string("dgp",DGPspec) * ".png"
sPathForOutput         = sFolderForOutput * sFileForOutput

savefig(sPathForOutput)
