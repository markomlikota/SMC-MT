# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Jan 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


Documentation for software folder for DSGE (RBC) illustration in
Mlikota & Schorfheide: "Sequential Monte Carlo With Model Tempering".

See "HowToSetupSMC-EstimationSoftware.txt" in "Documents/Notes"
for general notes on how best to set up SMC(-MT) estimation software.



# -------------------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------------------


In all of the following and in the codes in this folder (unless otherwise specified):

- "setup" refers to DGP and SMC specifications from DGPspec and SMCspec files

- "model" refers to Model1spec and Prior1spec files. Former contains fLL_Mi,
  latter fPriorDraw_Mi, fPriorLogEval_Mi, fIsDrawValid_Mi for the particular model Mi.

  e.g.  "Mi" = "L1" refers to first order linearized model,
        "Mi" = "L2" refers to second order linearized model,
        "Mi" = "L3" refers to third order linearized model,
        "Mi" = "VFI" refers to globally solved model (using value function iteration),

        "L11", "L21", "L31", "VFI1" are corresponding prior specifications where we estimate all parameters except ν
        "L12", "L22", "L32", "VFI2" estimate ν as well

- "estimation specification" refers to the way estimation is conducted
  (e.g. LTphiLast100, LTphiLast60, MTphiLast40_mL1_prL12_mc01)



# -------------------------------------------------------------------------------
# Julia Scripts
# -------------------------------------------------------------------------------


Many of the following scripts have several versions for different models, e.g. mL1, mL2, mL3, mVFI.


main_RBC_simulateData_mL1.jl :

simulates, stores and plots data for a given DGPspec under model "L1"

  - setting noMEs = true allows one to calibrate the StDs of measurement errors;
    it sets measurement errors to zero and returns the values one should take for ME StDs in the console
    in order to calibrate them to be a certain fraction of the overall observables' variances (specified by fracVarMEs)


main_RBC_OneRun.jl :

  carries out estimation of one model (as defined by Model1spec and Prior1spec)
  based on one SMC specification, DGP specification and using one estimation specification (M1-LT, or M0-LT, or M1-MT)

    - loads data based on DGPspec
    - loads model based on Model1spec and Prior1spec
    - in case of MT: loads M0 model based on Model0spec and Prior0spec
    - loads SMC tuning parameters based on SMCspec.
    - executes estimation (LT/MT)


main_RBC_RepeatedRuns.jl :

  carries out estimation of model "L2" using model "L1" as M0
  based on one SMC and DGP specification and using several estimation specifications and potentially many runs


main_RBC_outputAnalysisOneRun.jl :

  creates plots for output analysis of one DGP-SMC-Model1spec-Prior1spec and estimation specification (LT/MT)


main_RBC_outputAnalysisOneSpec.jl :

  creates plots for output analysis of one DGP-SMC-Model1spec-Prior1spec specification across many estimation specifications (LT/MT) and repeated runs


main_RBC_outputAnalysisMultipleSpecs.jl :

  for now, it plots only data under all (three) DGPs

  ... in future, it should also create plots that compare output under several DGP &/or SMC specs &/or estimation specs (LT/MT)


main_RBC_calibrateBSPF_mL2.jl :

  plots a histogram of LL evaluations under model "L2" for some DGPspec, using data simulated based on one of the models, given the number of particles to be used in BSPF


main_RBC_Checks_mL2.jl :

  performs some checks to see whether mL2 is coded correctly


main_RBC_ChecksBSPF.jl :

  compares LL evaluation under model L1 using the Kalman filter and the BSPF


main_RBC_LLevaltimes.jl :

  compares LL evaluation times under two models




# -------------------------------------------------------------------------------
# How to reproduce results (plots)
# -------------------------------------------------------------------------------


0) run "main_RBC_simulateData_mL2.jl" once using the following options:

  DGPspec     = "3"
  writeOutput = 1
  useModelspecInFileName = 2
  noMEs       = false
  fracVarMEs  = 0.05


1) run "main_RBC_RepeatedRuns.jl" once using the following options:

  DGPspec        = "3"
  SMCspec        = "03"
  vPhiLast       = [0.2, 0.4, 0.6, 0.8, 1.0]
  nRun           = 1


2) run "main_RBC_outputAnalysisOneSpec.jl" once using the following options:

  DGPspec        = "3"
  SMCspec        = "03"
  vPhiLast       = [0.2, 0.4, 0.6, 0.8, 1.0]
  nRun           = 1



# -------------------------------------------------------------------------------
# Output Storage
# -------------------------------------------------------------------------------



Folder "Data" contains

  - csv files of data simulated based on some DGPspec and Modelspec
  - corresponding plots of observables and latent states
  - plots of data under several DGPs and one Modelspec
  - plots of data under one DGPspec and several Modelspecs


Folder "Output": contains

  - a folder for each DGPspec-Model1spec-Prior1spec-SMC1spec combination
    - which contains a folder for each esitmation (LT/MT) specification (in case of MT: mentions Model0spec, Prior0spec, SMC0spec)
      - which contains the output of the SMC algorithm (csv files for first, last, all (and intermediate) stages)
      - and plots

  - a folder for various checks performed
