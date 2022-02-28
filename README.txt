# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Feb 2022, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


Documentation for overarching software folder for
Mlikota & Schorfheide: "Sequential Monte Carlo With Model Tempering".



# -------------------------------------------------------------------------------



This folder contains three sub-folders, each of which contain the software for
one of the three illustrations of SMC-MT that we consider in the paper;

  i) a simple example, using different univariate normal proposal distributions
      to sample from a univariate standard normal distribution

  ii) a VAR with stochastic volatility

  iii) a DSGE model: a simple RBC economy with capital adjustment costs

Each of these folders contains a separate readme, which shortly explains
the Julia scripts contained in the folder and contains instructions
on how to replicate results (plots) from the paper.



Also, this folder contains  a txt file "HowToSetupSMC-EstimationSoftware.txt",
which provides a best practice guide on how to set up software for estimating a
model using SMC-MT. It may be too general and will be revised in future.

Each of the above three folders contains a separate version of the SMC function
that allows for model tempering (although they are largely overlapping).
The best and most general function is in the folder for the RBC illustration.
