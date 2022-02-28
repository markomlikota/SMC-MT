# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Various helping functions.



# -------------------------------------------------------------------------------


function fGetCartProd(vSets,vIndices=0)

    # Computes Cartesian product of the sets in vSets with indices in vIndices

    # Inputs:  - vector of one-dimensional arrays containing the sets,
    #          - optional: vector of indices (one-dimensional Array) (default: all sets)

    # Output:  - matrix of points in resulting set (Cartesian product) across rows,
    #          - vector of tuples (all possible points in resulting space)

    if vIndices==0
        sets        = vSets
    else
        sets = [vSets[vIndices[i]] for i in 1:length(vIndices)]
    end

    tSpace      = unique(collect( Iterators.product(sets...) )) #tuple
    mSpace      = [tSpace[ind][i] for ind = 1:length(tSpace), i = 1:length(sets)] #matrix

    return mSpace, tSpace

end
