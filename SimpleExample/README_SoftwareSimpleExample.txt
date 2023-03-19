# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# March 2023, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


Documentation: software folder for "SimpleExample" illustration in
Mlikota & Schorfheide: "Sequential Monte Carlo With Model Tempering".


Section "How to reproduce results (plots)" below contains instructions on how to replicate our results. 
The other section is only meant to help understand the code/software-setup, but is not required for result-replication. 



# -------------------------------------------------------------------------------
# Julia Scripts
# -------------------------------------------------------------------------------


main_OneRun.jl :

Runs SMC once, given target and proposal distribution (& SMC settings).
These three are specified in folder "SpecFiles". Univariate and bivariate densities are supported.
Creates corresponding folder under "Output/" (based on target, proposal and SMC settings; e.g. tA1_pA1_smc1)
and stores SMC output (see fSetupSMC.jl) and plots.


main_MultipleProposalsAndRuns.jl :

Given target distribution (& SMC settings), runs SMC for many different proposal distributions, each many (nRun) times.
This code is only for the univariate case.
Creates corresponding folder under "Output/" (based on target and SMC settings; e.g. tA1_smc1),
in which the following are created:
  i)  a csv file with the parameters used for the different proposal distributions;
  ii) for each proposal distribution: a csv file with various estimation statistics and a csv file for the tempering schedule.
      These two files contain a row for each run.


main_MultipleProposalsAndRuns_d2.jl :

Same as the above, just for the bivariate case.


main_OutputAnalysis.jl :

Create plots based on output from "main_MultipleProposalsAndRuns.jl" (so for univariate case), stored in output folder created by latter.
Parts of it not written generically, so they work only for the exact options in "main_MultipleProposalsAndRuns.jl" that we used.


main_OutputAnalysis_d2.jl :

Same as the above, just for the bivariate case.



Note:

  - We used Targetspec and Propspec starting with "A" to indicate univariate dist. and "B" to indicate bivariate dist.

  - The functionality to create new folder "SetOfProposalsX+1" (i.e. append list) inside the output folder
    created by "main_MultipleProposalsAndRuns.jl" (e.g. "tA1_smc1") was created with the following intention.

    When one wants to add some more, different proposal parameters (means and/or stds and/or rhos, whatever),
    one should not have to re-run the estimation for the proposal parameter combinations already considered,
    but only for the newly added ones. For this the code would need to go to all the "SetOfProposalX" folders
    in the same specification path and determine which combinations have already been run, which are new.

    This is not implemented (yet). The code just runs all the parameter combinations under the new options
    considered and stores them in a new "SetOfProposals" folder. So for all practical purposes, just ignore this subfolder "SetOfProposalsX".



# -------------------------------------------------------------------------------
# How to reproduce results (plots)
# -------------------------------------------------------------------------------


1) run "main_oneRun.jl" once using the following options:

  Targetspec  = "A0"
  Propspec    = "A0"
  SMCspec     = "1"

2) run "main_MultipleProposalsAndRuns.jl" once using the following options:

  Targetspec          = "A0"
  SMCspec             = "1"
  nRun                = 100
  vMeans0             = -3:0.5:0
  vStds0              = 0.2:0.2:2
  setOfProposals      = 0


3) run "main_OutputAnalysis.jl" once using the following options:

  Targetspec      = "A0"
  SMCspec         = "1"
  setOfProposals  = 1
  vIndDM          = 3
