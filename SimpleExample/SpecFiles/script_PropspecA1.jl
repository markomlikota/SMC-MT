# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Specification of proposal distribution ( p(θ|Y,M0); (stage phiLast) posterior under M0)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


μg      = 0.0
σg      = 1.0


function fDrawProposal()

    return μg + randn() * σg

end


function fEvalLogProposal(θ)

    prob = pdf(Normal(μg,σg),θ)[1]

    return log(prob)

end


dProposal = Normal(μg,σg)
