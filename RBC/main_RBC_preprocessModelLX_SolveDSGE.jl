# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Jan 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Script to pre-process the model using package "SolveDSGE" (can be done once, before parameter values are defined, not for every θ anew), to be used then for Model "L2" (or "L3")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# LOAD PACKAGES: (NEED TO RUN ONLY ONCE, WHEN JULIA IS STARTED)
# -------------------------------------------------------------------------------


using SolveDSGE


# -------------------------------------------------------------------------------

# DEFINE PATH, INCLUDE EXTERNAL FUNCTIONS:
# -------------------------------------------------------------------------------


cd()
sMyPath            = string(pwd(),"/Dropbox/ResearchProjects/ModelTempering/Software/RBC/")
# sMyPath            = pwd() * "/"



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# --- MAIN BODY -----------------------------------------------------------------

# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


sPathModelspecSolveDSGE                = sMyPath * "SpecFiles/" * "script_ModelspecLX_SolveDSGE.txt"

process_model(sPathModelspecSolveDSGE)

sPathModelspecSolveDSGE_processed      = sMyPath * "SpecFiles/" * "script_ModelspecLX_SolveDSGE_processed.txt"
