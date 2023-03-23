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

