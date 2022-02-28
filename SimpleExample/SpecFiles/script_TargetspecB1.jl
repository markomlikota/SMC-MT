# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Specification of target distribution ( p(θ|Y,M1); posterior under M1)



# -------------------------------------------------------------------------------
# -------------------------------------------------------------------------------


μ1f      = 1.0
σ1f      = 0.4
μ2f      = 1.0
σ2f      = 0.8
ρf       = 0.4


covf     = ρf * σ1f * σ2f
mΣf      = [σ1f^2   covf    ;
            covf    σ2f^2   ]
vμf      = [μ1f,μ2f]


function fDrawTarget()

    return rand(MvNormal(vμf,mΣf))

end


function fEvalLogTarget(θ)

    prob = pdf(MvNormal(vμf,mΣf),θ)[1]

    return log(prob)

end


dTarget = MvNormal(vμf,mΣf)
