# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Specification of target distribution ( p(θ|Y,M1); posterior under M1)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


μf      = 1.0
σf      = 0.2


function fDrawTarget()

    return μf + randn() * σf

end


function fEvalLogTarget(θ)

    prob = pdf(Normal(μf,σf),θ)[1]

    return log(prob)

end


dTarget = Normal(μf,σf)
