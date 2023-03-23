# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Various helping functions related to (implementation of) SMC algorithm.



# -------------------------------------------------------------------------------






function fGetParaPointers(vGroupLengths)

    # Creates X x 2 matrix of parameter pointers that facilitate going from a vector of parameters vθ to the X groups it is composed of; each row indicates the first and last index that has to be used in order to go from vθ to the respective group:
    # e.g. Φ = reshape(θ[mParaPointers[1,1]:mParaPointers[1,2]],k,n)

      vNParams              = Int.(vGroupLengths)
      mParaPointers         = Int.([ones(length(vNParams)) vNParams])
      [mParaPointers[nn,:]  .+= mParaPointers[nn-1,2] for nn = 2:length(vNParams)]

      return mParaPointers

end




function fMyGamma(mean,std)

    # Computes the two parameters for Gamma distribution, given desired mean and std.

    var = std^2
    θ   = var/mean
    α   = mean/θ

    return α, θ

end

function fMyInvGamma(v,s)

    # Goes from specification in HS(15),AS(07) to that of Julia & Wikipedia.
    # Note: need to divide draw by 100.

    α   = v
    θ   = v*s^2 / 2

    return α, θ

end

function fMyInvGamma2(p1,p2)

    # Goes from specification in AMSV to that of Julia & Wikipedia.
    # Note: need to take square root of draw.

    a = p2/2
    b = p1^2 * a

    return a,b

end

function fMyBeta(mean,std)

    # Computes the two parameters for Gamma distribution, given desired mean and std.

    var = std^2

    α   = (1-mean)*mean^2 / var + mean
    β   = α * (1-mean)/mean

    return α, β

end

# myChol(x) = cholesky(Matrix(Hermitian(x))).L
