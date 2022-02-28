# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# DGP specification file:

# "DGP 2"; T=100, high ξ, lower ρ (than "DGP 1")



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------

# Sample size:

T     = 100


# True parameter-matrices:

Φ10   = [0.6 0.3;
         0.0 0.4]
Φc0   = [0.0;0.0]
Φ0    = [Φ10 Φc0]'

Φε0   = [1.0 0.0;
         0.7 1.0] # because fSimulateVARSVp1 uses Φε0, not Σ0
Σ0    = Φε0 * Φε0'

vd0   = [1, 1]
DD0   = diagm(vd0)

vρ0   = [0.20, 0.60]
vξ0   = [0.80, 0.90]

n     = 2
p     = 1
k     = n*p +1
