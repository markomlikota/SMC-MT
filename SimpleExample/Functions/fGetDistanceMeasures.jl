# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Function to get distance/overlap measures of proposal and target.



# -------------------------------------------------------------------------------


function solvequadratic(a, b, c)
    d = sqrt(b^2 - 4a*c)
    (-b - d) / 2a, (-b + d) / 2a
end



function fGetDistanceMeasures(vμ1,mΣ1,vμ0,mΣ0)


    k           = length(vμ1)


    # Area of product of two normals:

    P1          = inv(mΣ1)
    P0          = inv(mΣ0)

    μbar        = inv(P0+P1)* (P0* vμ0 + P1* vμ1)
    Oterm1      = (2*π)^(-1/2) * ( det(P0)*det(P1)/det(P0+P1) )^(1/2)
    Oterm2      = exp( -1/2 * ( vμ0'* P0* vμ0 + vμ1'* P1* vμ1 - μbar'* (P0+P1)* μbar ) )

    O           = Oterm1 * Oterm2
    AREAPROD    = O

    # Kullback-Leibler distance:

    KLterm1     = logdet(P1) - logdet(P0)
    KLterm2     = tr( P1*(inv(P0) + vμ0* vμ0') - 2* P1* vμ1* vμ0' + P1* vμ1* vμ1' - I)
    KL          = 1/2 * ( KLterm1 - KLterm2 )

    # Variance (under proposal) of ratio of target-pdf to proposal-pdf (from Importance Sampling):

    # Analytical expression obtained only for proposal variance at least as high as half the target variance (and analogously for bivariate...):
    mH          = 2* inv(mΣ1) - inv(mΣ0)
    mK          = 2* inv(mΣ1) * vμ1 - inv(mΣ0) * vμ0
    if sum( eigvals(mH) .> 0 ) == k
        ISterm1     = sqrt( det(mΣ0) )/det(mΣ1) * sqrt( det(inv(mH)) )
        ISterm2     = exp( -1/2 * ( vμ1' * 2*inv(mΣ1)*vμ1 - vμ0' * inv(mΣ0)*vμ0 - mK'*inv(mH)*mK ) )
        IS          = ISterm1 * ISterm2 - 1
        IS          = log(1+IS)
    else
        IS          = NaN
    end

    # Area of minimum of two normals:

    if k == 1

        if vμ0 < vμ1
            μa = vμ0
            μb = vμ1
            σa = sqrt(mΣ0)
            σb = sqrt(mΣ1)
        else
            μa = vμ1
            μb = vμ0
            σa = sqrt(mΣ1)
            σb = sqrt(mΣ0)
        end

        termA = 1/2* (1/σb^2 - 1/σa^2)
        termB = (μa/σa^2 - μb/σb^2)
        termC = log(σb/σa) + 1/2* ( μb^2/σb^2 - μa^2/σa^2 )

        if termA == 0
            vc = repeat([-termC/termB],2,1)
        else
            vc     = solvequadratic(termA,termB,termC)
        end

        c1 = minimum(vc)
        c2 = maximum(vc)

        if c1 == c2
            AREAMIN = normcdf(μb,σb,c1) + 1-normcdf(μa,σa,c1)
        else
            if σa < σb
                AREAMIN = 1-normcdf(μa,σa,c2) + normcdf(μa,σa,c1) + normcdf(μb,σb,c2) -normcdf(μb,σb,c1)
            else
                AREAMIN = 1-normcdf(μb,σb,c2) + normcdf(μb,σb,c1) + normcdf(μa,σa,c2) -normcdf(μa,σa,c1)
            end
        end

        if μa == μb && σa == σb
            AREAMIN = 1
        end

    else

        AREAMIN = NaN

    end


    return [AREAPROD, KL, 1-AREAMIN, IS]

end
