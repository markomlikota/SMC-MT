# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Function to simulate VAR(1) with Stochastic Volatility. Only works for p=1 !



# -------------------------------------------------------------------------------



function fSimulateVARSVp1(Φ,Φε,vd,vρ,vξ,T,vy0=0,vdt0=0)


      # Extract all objects needed:

      Φ1          = Φ'[:,1:end-1]
      Φc          = Φ'[:,end]

      n           = length(Φc)


      # Initialize y, d:

      vy0 == 0 ? vy0 = inv(I-Φ1)*Φc : nothing
      vdt0 == 0 ? vdt0 = ones(n,1) : nothing


      # Simulate d_t:

      mη                = randn(n,T)

      mLogdt            = zeros(n,T+1)
      mLogdt[:,1]       = log.(vdt0)
      [mLogdt[:,t+1]    = vρ .* mLogdt[:,t] + vξ .* mη[:,t] for t=1:T]

      mLogd             = mLogdt .+ log.(vd)
      md                = exp.(mLogd)


      # Create objecsts to store y_t, Σ_t, ε_t, and initialize them:

      my                = zeros(n,T+1)
      my[:,1]           = vy0

      vmΣ               = Matrix{Float64}[zeros(n,n) for t in 1:(T+1)]
      vmΣ[1]            = Φε * diagm(md[:,1]) * Φε'

      mε                = zeros(n,T+1)


      # Simulate for t = 1,2,...,T:

      for t = 2:(T+1)
            mDt         = diagm(md[:,t])
            vmΣ[t]      = Φε * mDt * Φε'
            ε           = rand(MvNormal(zeros(n),mDt))
            my[:,t]     = Φ1 * my[:,t-1] + Φc + Φε * ε
            mε[:,t]     = ε
      end


      return my', vmΣ, md

end
