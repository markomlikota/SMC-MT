# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Functions to work with VARs.



# -------------------------------------------------------------------------------



function fDataVARtoReg(mData,p)

    # Creates matrices Y (T x n) and X (T x k) from ((T+p) x n) matrix of data.
    # This allows to write VAR(p) as Y = X Φ + U.

    n           = size(mData)[2]

    Y           = mData[p+1:end,:]

    T           = size(Y)[1]
    k           = n*p+1
    X           = zeros(T,k)
    [X[t,:]     = [vec(mData[ reverse( ((t-1) .+(1:p)) ) ,:]'); 1]' for t=1:T]

    return Y, X

end


function fCompanionForm(Φ, Σ)

    # Computes matrices F0, F1, Fε, M for companion form (to write VAR(p) as VAR(1)), given Φ and Σ (as defined in ECON 722 slides).

    k, n            = size(Φ)
    p               = Int((k-1)/n)

    F1              = zeros(n*p,n*p)
    [F1[1:n,(pp-1)*n.+(1:n)]    = Φ[(pp-1)*n.+(1:n),1:n]' for pp=1:p]
    F1[n+1:n*p,1:n*p-n]         = diagm(vec(ones(n*p-n,1)))

    F0              = zeros(n*p,1)
    F0[1:n]         = Φ[end,:]

    Fε              = zeros(n*p,n*p)
    Fε[1:n,1:n]     = cholesky(Σ).L

    M               = [diagm(vec(ones(n,1))) zeros(n,n*p-n)]

    return F0, F1, Fε, M

end


function fPosteriorMDDunderVARMNIW(Y,X,μ_pr,P_pr,v_pr,S_pr,ϕ_Nϕ=1.0)

      # Computes posterior distribution parameters and logMDD for a VAR with MNIW prior

      # Inputs: data, prior MNIW parameters, "early stopping parameter" ϕ_Nϕ

      T,n  = size(Y)
      k    = size(X,2)

      P_po = P_pr  + ϕ_Nϕ * X' * X
      μ_po = inv(P_po) * ( P_pr * μ_pr  + ϕ_Nϕ * X'*Y )


      S_po = S_pr  + μ_pr' * P_pr * μ_pr  + ϕ_Nϕ * Y'*Y  - μ_po' * P_po * μ_po
      v_po = v_pr  + ϕ_Nϕ* T

      logMDD = -T*ϕ_Nϕ*n/2 * log(π) + n/2* (-logdet(P_po) + logdet(P_pr) ) - v_po/2 * logdet(S_po) + v_pr/2 * logdet(S_pr) +  logmvgamma(n,v_po/2) - logmvgamma(n,v_pr/2)

      return μ_po, P_po, v_po, S_po, logMDD

end
