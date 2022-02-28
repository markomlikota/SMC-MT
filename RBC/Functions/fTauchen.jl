# Marko Mlikota, University of Pennsylvania, mlikota@sas.upenn.edu
# July 2020
# -------------------------------------------------------------------------------



function fTauchen(na,abar,ρ,σ,l=-Inf,u=Inf,own_a_bounds=[0 0])

    # Applies Tauchen(1986)-Method to construct grid for exogenous state variable a and its transition matrix.
    # Process is         a' = (1-ρ)abar + ρa + σ ε

    # Inputs:   - number of gridpoints "na" (odd)
    #           - mean value "abar"
    #           - persistence parameter "ρ"
    #           - standard deviation of "overall" error term (σε) "σ"
    #           - lower and upper bound for error term (ε) (if it follows a truncated normal)
    #           - 1x2 vector of own bounds for distance of mean to grid-ends expressed in units of variance (optional; default bounds: +/- m = (na-1)/2)

    # Output:   - vector of values
    #           - transition matrix mP



    if own_a_bounds != [0 0]
        va = range(own_a_bounds[1],own_a_bounds[2],length=na)
        va = collect(va)
        va = abar .+ va*σ/sqrt(1-ρ^2)
    else
        m = (na-1)/2;
        a_1 = abar - m*σ/sqrt(1-ρ^2);
        a_na = abar + m*σ/sqrt(1-ρ^2);
        va = range(a_1,a_na, length=na); #list of na linearly spaced points from a_1 to a_5
        va = collect(va);
    end
    Δ_a = va[2] - va[1];

    mP = zeros(na,na);
    for i in 1:na, j in 1:na
        μ = (1-ρ)*abar + ρ*va[i];
        if j == 1

            # mP[i,j] = cdf(TruncatedNormal(μ, σ, μ+σ*l, μ+σ*u), va[j] + Δ_a/2) - 0;
            mP[i,j] = cdf(truncated(Normal(μ, σ), μ+σ*l, μ+σ*u), va[j] + Δ_a/2) - 0
        elseif j == na
            # mP[i,j] = 1 - cdf(TruncatedNormal(μ, σ, μ+σ*l, μ+σ*u), va[j] - Δ_a/2);
            mP[i,j] = 1 - cdf(truncated(Normal(μ, σ), μ+σ*l, μ+σ*u), va[j] - Δ_a/2)
        else
            # mP[i,j] = cdf(TruncatedNormal(μ, σ, μ+σ*l, μ+σ*u), va[j] + Δ_a/2) - cdf(TruncatedNormal(μ, σ, μ+σ*l, μ+σ*u), va[j] - Δ_a/2);
            mP[i,j] = cdf(truncated(Normal(μ, σ), μ+σ*l, μ+σ*u), va[j] + Δ_a/2) - cdf(truncated(Normal(μ, σ), μ+σ*l, μ+σ*u), va[j] - Δ_a/2);
        end
    end

    return va, mP
end
