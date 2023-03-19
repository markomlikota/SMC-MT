# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# March 2023, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


Documentation: software folder for
Mlikota & Schorfheide: "Sequential Monte Carlo With Model Tempering".



# -------------------------------------------------------------------------------



This folder contains four sub-folders, each of which contain the software for
one of the four illustrations of SMC-MT that we consider in the paper;

  i) a simple example, using different univariate normal proposal distributions
      to sample from a univariate standard normal distribution

  ii) a VAR with stochastic volatility

  iii) a DSGE model: a simple RBC economy with capital adjustment costs

  iv) a larger DSGE model: New Keynesian DSGE model from Aruoba, Bocola & Schorfheide (2017)

Each of these folders contains a separate readme, which shortly explains
the Julia scripts contained in the folder and contains instructions
on how to replicate results (plots) from the paper.



Also, this folder contains  a txt file "HowToSetupSMC-EstimationSoftware.txt",
which provides a best practice guide on how to set up software for estimating a
model using SMC-MT.

Each of the above three folders contains a separate version of the SMC function
that allows for model tempering (although they are largely overlapping).
The best and most general function is in the folder for the RBC illustration.
