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

Sidenote: Each of these folders contains a separate version of the SMC function
(although they are largely overlapping). 
The best and most general function is in the folder for the RBC illustration.



Also, there is a fifth folder "Template" containing the template-script "SMCMT.jl", 
which can be easily adjusted to apply SMC-MT for one's own application, 
with accompanying instructions and comments. 
The script "SMCMT_example.jl" uses this template to estimate a very simple model.


Depending on one's computer, one might need to adjust the functions writing, reading and appending cdv-files.
See comment and code in fFolderFileManagment.jl, which exists in each of application-folders.



