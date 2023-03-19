# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# March 2023, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


Documentation for software folder for second DSGE illustration (ABS 2017 model) in
Mlikota & Schorfheide: "Sequential Monte Carlo With Model Tempering".

Section "How to reproduce results (plots)" below contains instructions on how to replicate our results. 
The other sections are only meant to help understand the code/software-setup, but are no required for result-replication. 



# -------------------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------------------


In all of the following and in the codes in this folder (unless otherwise specified):

- "setup" refers to DGP and SMC specifications from DGPspec and SMCspec files

- "model" refers to Model1spec and Prior1spec files. Former contains fLL_Mi,
  latter fPriorDraw_Mi, fPriorLogEval_Mi, fIsDrawValid_Mi for the particular model Mi.

  e.g.  "Mi" = "L1" refers to first order linearized model,
        "Mi" = "L2" refers to second order linearized model
	
	each has a corresponding prior specification file

- "estimation specification" refers to the way estimation is conducted
  (e.g. LTphiLast100, MTphiLast100_mL1_prL1_mc1)



# -------------------------------------------------------------------------------
# Julia Scripts
# -------------------------------------------------------------------------------


main_ABS_simulateData_mL1.jl :
  
  simulates, stores and plots data for a given DGPspec under model "L1", 
  analogous file for model "L2"

  - setting noMEs = true allows one to calibrate the StDs of measurement errors;
    it sets measurement errors to zero and returns the values one should take for ME StDs in the console
    in order to calibrate them to be a certain fraction of the overall observables' variances (specified by fracVarMEs)


main_ABS_OneRun.jl :

  carries out estimation of one model (as defined by Model1spec and Prior1spec)
  based on one SMC specification, DGP specification and using one estimation specification (M1-LT, or M0-LT, or M1-MT)

    - loads data based on DGPspec
    - loads model based on Model1spec and Prior1spec
    - in case of MT: loads M0 model based on Model0spec and Prior0spec
    - loads SMC tuning parameters based on SMCspec.
    - executes estimation (LT/MT)


main_ABS_outputAnalysisOneRun.jl :

  creates plots for output analysis of one DGP-SMC-Model1spec-Prior1spec and estimation specification (LT/MT)


main_ABS_calibrateBSPF_mL2.jl :

  plots a histogram of LL evaluations under model "L2" for some DGPspec, using data simulated based on one of the models, given the number of particles to be used in BSPF


main_ABS_Checks_mL2.jl :

  performs some checks to see whether mL2 is coded correctly


main_ABS_preprocessModelLX_SolveDSGE.jl :
 
  Script to pre-process model written as .txt file so that package SolveDSGE can perform second- (and first-)order linearizations



# -------------------------------------------------------------------------------
# How to reproduce results (plots)
# -------------------------------------------------------------------------------


0) run "main_ABS_simulateData_mL2.jl" once using the following options:

  DGPspec     = "1"
  writeOutput = 1
  useModelspecInFileName = 2
  noMEs       = false
  fracVarMEs  = 0.05


1) run "main_ABS_OneRun.jl" once using the following options:

	DGPspec         = "1"     
	Model1spec      = "L2"    
	Prior1spec      = "L2"    
	SMC1spec        = "1"    

	sTemperingType  = "LT"         
	ϕ_Nϕ            = 1.0

	          
2) run "main_ABS_OneRun.jl" once using the following options:

	DGPspec         = "1"     
	Model1spec      = "L1"    
	Prior1spec      = "L1"    
	SMC1spec        = "1"    

	sTemperingType  = "LT"         
	ϕ_Nϕ            = 1.0


3) run "main_ABS_OneRun.jl" once using the following options:

	DGPspec         = "1"    
	Model1spec      = "L2"    
	Prior1spec      = "L2"   
	SMC1spec        = "1"    

	sTemperingType  = "MT"       
	ϕ_Nϕ            = 1.0       


	Model0spec      = "L1"        
	Prior0spec      = "L1"        
	ϕ_Nϕ_M0         = 1.0    
	SMC0spec        = "1"  


Note: 
	- can interrupt run and continue later. For this, when taking up the run again, set "continueInterruptedRun = 1". 
	- preprocessing of model required for SolveDSGE package already performed. 
	  This involves simply running the above mentioned script "main_ABS_preprocessModelLX_SolveDSGE.jl" as-is.



# -------------------------------------------------------------------------------
# Output Storage
# -------------------------------------------------------------------------------



Folder "Data" contains

  - csv files of data simulated based on some DGPspec and Modelspec
  - corresponding plots of observables and latent states
  - plots of data under one DGPspec and both Modelspecs


Folder "Output": contains

  - folder "dgp1_mL1_prL1_mc1/LTphiLast100" with results for SMC-LT estimation of model "L1" (output of the SMC algorithm (and plots, if output-analysis-script is run))

  - folder "dgp1_mL2_prL2_mc1/LTphiLast100" with results for SMC-LT estimation of model "L2" (output of the SMC algorithm (and plots, if output-analysis-script is run))

  - folder "dgp1_mL2_prL2_mc1/MTphiLast100_mL1_prL1_mc1" with results for SMC-MT estimation of model "L2" (output of the SMC algorithm (and plots, if output-analysis-script is run))

	Note: to get the runtime reduction reported in the paper, add to the time of the SMC-MT run of L2 the time of the preliminary SMC-LT run of L1.

  - a folder for various checks performed
