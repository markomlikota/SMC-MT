# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Function to draw from multivariate normal distribution when Σ is not positive definite



# -------------------------------------------------------------------------------


function fMyMNDraw(vμ,mΣ)

      n = length(vμ)
      F = svd(mΣ)

      return vμ + F.U * Diagonal(sqrt.(F.S)) * rand(MvNormal(zeros(n),I))

end
