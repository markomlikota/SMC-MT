# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Specification of target distribution ( p(θ|Y,M1); posterior under M1)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


μf      = 0.0
σf      = 1.0


function fDrawTarget()

    return μf + randn() * σf

end


function fEvalLogTarget(θ)

    # prob = pdf(Normal(μf,σf),θ)[1]
    prob = logpdf(Normal(μf,σf),θ)[1]

    return prob

end


dTarget = Normal(μf,σf)
