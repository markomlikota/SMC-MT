# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# March 2023, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


Documentation: software folder for VAR illustration in
Mlikota & Schorfheide: "Sequential Monte Carlo With Model Tempering".

Section "How to reproduce results (plots)" below contains instructions on how to replicate our results. 
The other sections are only meant to help understand the code/software-setup, but are no required for result-replication. 



# -------------------------------------------------------------------------------
# Preliminary
# -------------------------------------------------------------------------------


In all of the following and in the codes in this folder (unless otherwise specified):

- "setup" refers to DGP and SMC specifications from DGPspec and SMCspec files

        "DGPspec1", "DGPspec2" and "DGPspec3" refer to DGP A, B and C, respectively
        "SMCspec1" and "SMCspec2" refer to N=500 and N=1000 particles, respectively

- "model" refers to Model1spec and Prior1spec files. Former contains fLL_Mi,
  latter fPriorDraw_Mi, fPriorLogEval_Mi, fIsDrawValid_Mi for this particular model Mi.

  e.g.  "Mi" = "HM" refers to homoskedastic VAR(1) with intercept,
        "Mi" = "SV" refers to VAR(1) with Stochastic Volatility,
        and "HM1" and "SV1" are corresponding prior specifications based on Minnesota prior for Φ, Σ

- "estimation specification" refers to the way estimation is conducted
  (e.g. LTphiLast100, LTphiLast60, MTphiLast40_mHM_prHM1)



# -------------------------------------------------------------------------------
# Julia Scripts
# -------------------------------------------------------------------------------


main_VARSV_simulateData.jl :

simulates, stores and plots data for a given DGPspec


main_VARSV_OneRun.jl :

  carries out estimation of one model (as defined by Model1spec and Prior1spec)
  based on one setup (DGP & SMC specification) and using one estimation specification (M1-LT, or M0-LT, or M1-MT)

    - loads data based on DGPspec
    - loads model based on Model1spec and Prior1spec
    - in case of MT: loads M0 model based on Model0spec and Prior0spec
    - loads SMC tuning parameters based on SMCspec.
      Note: ϕ_Nϕ specified separately on top, in preamble "Options", not in SMCspec.

    - executes estimation:
      - under LT: takes vParLabs, fLL, fLLtilde, fPriorDraw, fPriorLogEval, fIsDrawValid
        from Model1spec and Prior1spec
      - under MT: first defines all of the above with suffixes M1 and M0 based on
        (Model1spec and Prior1spec) and (Model0spec and Prior0spec) respectively.
        Then uses these to define the above without suffixes based on how
        the common parameter vector θ is assembled from ϑM1 and ϑM0.


main_VARSV_RepeatedRuns.jl :

  carries out estimation of one model, based on one setup for several estimation specifications and all this for many runs.

    - loads data based on DGPspec
    - loads SMC tuning parameters based on SMCspec

    - executes estimation of M1-LT, M0-LT and M1-MT, respectively, many (nRun) times,
      latter two each time for several ϕ_Nϕ_M0's (vector specified in preamble)


main_VARSV_RepeatedRuns_v2.jl :

  an improved version of the above, where the models to be estimated can be specified explicitly.
  In the above, the run folders are named "runX", and here there is the option to choose a name.
  Thus, in the runs of the actual model SV under different DGPs and SMCspecs, we have the folders "runX".
  For the runs of the models SVxX with artificially increased runtimes, we chose the name "xrunX".


main_VARSV_outputAnalysisOneSpec.jl :

  creates plots for output analysis of one DGP & SMC spec over all estimation specifications and many runs

  Some parts are not generically coded; they work only for our considered Modelspecs "HM" and "SV" (and Priorspecs)


main_VARSV_outputAnalysisMultipleSpecs.jl :

  creates plots for output analysis that combine several DGP &/or SMC specs


  main_VARSV_outputAnalysisMultipleSpecsLLTimes.jl :

    analogous to the one above, just for the models which are equal to SV, but with artificially incrased runtimes


main_VARSV_LLevaltimes.jl :

  compares LL evaluation times under different models


main_testSMC.jl :

  tests SMC algorithm using simple example (analytical solution available)

  Note:
    - running main_VARSV_OneRun.jl for Model1spec = "VAR" or "HM" and corresponding Prior1specs "VAR1" and "HM1"
      can also serve as test; analytical solution is computed at bottom of that script.



# -------------------------------------------------------------------------------
# How to reproduce results (plots)
# -------------------------------------------------------------------------------


0) run "main_VARSV_simulateData.jl" three times with the following options, respecively:

  i)    DGPspec   = "1"
  ii)   DGPspec   = "2"
  iii)  DGPspec   = "3"


1) run "main_VARSV_RepeatedRuns.jl" four times with the following options, respectively:

  i)    DGPspec = "1", SMCspec = "1"
  ii)   DGPspec = "2", SMCspec = "1"
  iii)  DGPspec = "3", SMCspec = "1"
  iv)   DGPspec = "1", SMCspec = "2"

  For each, set
    vPhiLast        = [0.2, 0.4, 0.6, 0.8, 1.0]
    nRun            = 200

  Note:
    - on a 8-core, 16gb Macbook Pro 2016, one run under each of the above takes
      about 30min, 45min, 80min, 60min, respectively.
    - one can also get to 200 runs by running the code several times,
      each time setting nRun to some smaller number and running the script anew.
      Code will automatically continue adding new run folders where it stopped last time.


2) run "main_VARSV_outputAnalysisOneSpec.jl" four times using the same options as in 1) above.


3) run "main_VARSV_outputAnalysisMultipleSpecs.jl" once using the following options:

      vDGPspecs         = ["1","2","3","1"]
      vSMCspecs         = ["1","1","1","2"]
      vPhiLast          = [0.2, 0.4, 0.6, 0.8, 1.0]
      nRun              = 200


4) run "main_VARSV_RepeatedRuns_v2.jl" four times with the following options, respectively:

  i)    Model1spec = "SV",    Prior1spec = "SV1"
  ii)   Model1spec = "SVx2",  Prior1spec = "SVx21"
  iii)  Model1spec = "SVx4",  Prior1spec = "SVx41"

  For each, let the other options be:

    DGPspec         = "1"
    SMCspec         = "1"
    Model0spec      = "HM"
    Prior0spec      = "HM1"
    vPhiLast        = [0.2, 0.4, 0.6, 0.8, 1.0]
    nRun            = 1
    sRunFolderName  = "xrun"
    seed            = 2


5) run "main_VARSV_outputAnalysisMultipleSpecsLLTimes.jl" once using the following options:

  vModel1specs            = ["SV","SVx2","SVx4","SVx8"]
  vPrior1specs            = ["SV1","SVx21","SVx41","SVx81"]
  DGPspec                 = "1"
  SMCspec                 = "1"
  Model0spec              = "HM"
  Prior0spec              = "HM1"
  vPhiLast                = [0.2, 0.4, 0.6, 0.8, 1.0]
  nRun                    = 1
  sRunFolderName          = "xrun"



# -------------------------------------------------------------------------------
# Output Storage
# -------------------------------------------------------------------------------


Folder "Output": contains

- a folder for each DGP-SMC-Model-Prior specification
- all plots that combine several DGP &/or SMCs


Each DGP-SMC-Model-Prior specification folder contains

- a folder "OneRun" created by main_VARSV_OneRun.jl (this was preliminary testing, running a specification once)
- many "runX" folders and/or "xrunX" folders
- all plots that relate only to this DGP-SMC


In case of model "HM", there are folders named "LTphiLast60" (and other values for phiLast, including "LTphiLast100", full LT estimation).
In case of model "SV", there are is a folder named "LTphiLast100" and folders named "MTphiLast40_mHM_prHM1" for different values of phiLast.

Each of these folders, in turn, contain the actual SMC output (see preamble of fSetupSMC.jl).
