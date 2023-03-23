# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Function to get proposal objects.



# -------------------------------------------------------------------------------




function fGetProposalObjects(μg,σg)


    function fDrawProposal()

        return μg + randn() * σg

    end


    function fEvalLogProposal(θ)

        prob = logpdf(Normal(μg,σg),θ)[1]

        return prob

    end


    dProposal = Normal(μg,σg)

    return fDrawProposal, fEvalLogProposal, dProposal

end



function fGetProposalObjects_d2(μ1g,σ1g,μ2g,σ2g,ρg)

    covg     = ρg * σ1g * σ2g
    mΣg      = [σ1g^2   covg    ;
                covg    σ2g^2   ]
    vμg      = [μ1g,μ2g]


    function fDrawProposal()

        return rand(MvNormal(vμg,mΣg))

    end


    function fEvalLogProposal(θ)

        prob = logpdf(MvNormal(vμg,mΣg),θ)[1]

        return prob

    end


    dProposal = MvNormal(vμg,mΣg)

    return fDrawProposal, fEvalLogProposal, dProposal

end
