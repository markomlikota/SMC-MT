# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Specification of proposal distribution ( p(θ|Y,M0); (stage phiLast) posterior under M0)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


μ1g      = 0.0
σ1g      = 1.0
μ2g      = 0.0
σ2g      = 1.0
ρg       = 0.0


covg     = ρg * σ1g * σ2g
mΣg      = [σ1g^2   covg    ;
            covg    σ2g^2   ]
vμg      = [μ1g,μ2g]


function fDrawProposal()

    return rand(MvNormal(vμg,mΣg))

end


function fEvalLogProposal(θ)

    prob = pdf(MvNormal(vμg,mΣg),θ)[1]

    return log(prob)

end


dProposal = MvNormal(vμg,mΣg)
