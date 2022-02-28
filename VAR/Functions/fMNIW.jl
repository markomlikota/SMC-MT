# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Functions to draw from and evaluate MNIW.



# -------------------------------------------------------------------------------


function fOneDrawMNIW(μ,P,S,v)

    Σ     = rand( InverseWishart(v,Matrix(Hermitian(S))) )
    Φ     = rand( MatrixNormal(μ, Matrix(Hermitian(inv(P))), Matrix(Hermitian(Σ)) ))

    return Σ, Φ

end

function fOneDrawMNIWstationary(μ,P,S,v)

    # Draws one Phi and Sigma, given (posterior) parameters (MNIW()).

    while true

        Σ     = rand( InverseWishart(v,Matrix(Hermitian(S))) )
        Φ     = rand( MatrixNormal(μ, Matrix(Hermitian(inv(P))), Matrix(Hermitian(Σ)) ))

        maximum(abs.(eigvals(fCompanionForm(Φ,Σ)[2]))) < 1 ? (return Σ, phi) : nothing

    end

end

function fEvalMNIW(Φ,Σ,μ,P,S,v)

    prprΣ   = pdf(InverseWishart(v,S),Σ)
    prprΦ   = pdf( MatrixNormal(μ, Matrix(Hermitian(inv(P))), Matrix(Hermitian(Σ)) ), Φ )

    return prprΣ * prprΦ

end


#
# function fDrawMNIW(μBar,pBar,sBar,vBar,nDraws)
#
#     # Draws Phi and Sigma, given (posterior) parameters (MNIW()).
#
#     @eval @everywhere μBar = $μBar
#     @eval @everywhere pBar = $pBar
#     @eval @everywhere sBar = $sBar
#     @eval @everywhere vBar = $vBar
#
#     @everywhere k       = size(μBar)[1]
#     @everywhere n       = size(sBar)[2]
#     aSig    = SharedArray{Float64}(nDraws,n,n)
#     aPhi    = SharedArray{Float64}(nDraws,k,n)
#
#     @sync @distributed for dd = 1:nDraws
#
#         sig, phi = fOneDrawMNIWstationary(μBar,pBar,sBar,vBar)
#
#         aSig[dd,:,:]    = sig
#         aPhi[dd,:,:]    = phi
#     end
#
#     return aSig, aPhi
#
# end
