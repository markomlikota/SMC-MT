# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Functions for prior-posterior plot and bridge distribution-evolution plot



# -------------------------------------------------------------------------------


function fFindPlotBounds(boundsType,vx)

    # Finds plot bounds from vector vx obtained via Kernel Density Estimate function "kde" given desired boundsType

    # boundsType is either 0 (whole real line), 1 (positive real line), or 2 ((0,1) interval)

    if boundsType == 0
        st = 1; en = length(vx)
    else
        st = findfirst(vx .>= 0)
        boundsType == 1 ? en = length(vx) : en = findlast(vx .<= 1)
    end

    return st, en

end

function fPrPoPlot(pInd,aParticles,mWeights,boundsType)

    # Prior-Posterior Plot

    # Inputs:   - pInd is index of parameter
    #           - (N x np x Nϕ) array of particles
    #           - (N x Nϕ) matrix of weights
    #           - boundsType is either 0 (whole real line), 1 (positive real line), or 2 ((0,1) interval)

    k       = kde(aParticles[:,pInd,end], weights = Weights(mWeights[:,end]))
    st, en  = fFindPlotBounds(boundsType,k.x)
    ppp     = fMyPlot(k.x[st:en],k.density[st:en],vParLabs[pInd])

    k       = kde(aParticles[:,pInd,1], weights = Weights(mWeights[:,end]))
    st, en  = fFindPlotBounds(boundsType,k.x)
    plot!(k.x[st:en],k.density[st:en],label="",line=(:black,0.9,2,:dot),yticks=false)

    return ppp

end

# pInd,yLab,aParticles,mWeights,nLines,boundsType = pInd,yLab,aAllParticlesM0[:,:,2:6],mWeights,nLines,boundsType
function fPlotDistEvol(pInd,yLab,aParticles,mWeights,nLines,boundsType)

    # Plots evolution of bridge distributions

    # Inputs:   - pInd is index of parameter
    #           - (N x np x Nϕ) array of particles
    #           - (N x Nϕ) matrix of weights
    #           - nLines is number of bridge distributions to plot
    #           - boundsType is either 0 (whole real line), 1 (positive real line), or 2 ((0,1) interval)



    Nϕ          = size(aParticles)[3]

    vn          = Int.( round.( collect(range(1,Nϕ,length = nLines)) ))

    transpStart = 0.5
    transpEnd   = 1
    vTransp     = collect(range(transpStart,transpEnd,length=nLines))

    k           = kde(aParticles[:,pInd,1], weights = Weights(mWeights[:,end]))
    st, en      = fFindPlotBounds(boundsType,k.x)


    vx          = k.x[st:en] #range(-0.2,1.0,length=2048) #
    vy          = k.density[st:en] #range(0.0,2.0,length=2048) #
    ppp         = plot(ones(en-st+1),vx,vy, camera=(60,30), line=(:black,transpStart,2,:line),label="",xlabel="n",ylabel=yLab,ylims=(k.x[st],k.x[en]),
    xtickfont=font(16),ytickfont=font(16),guidefont=font(16),legendfontsize=16 )
    xflip!(true)


    for nn in 2:nLines

        n       = vn[nn]

        k       = kde(aParticles[:,pInd,n], weights = Weights(mWeights[:,end]))
        st, en  = fFindPlotBounds(boundsType,k.x)

        ppp     = plot!(ones(en-st+1)*n,k.x[st:en], k.density[st:en],camera=(60,30),line=(:black,vTransp[nn],2,:line),label="")

    end


    return ppp

end


function fPlotDistEvol2(pInd,yLab,aParticles,mWeights,nLines,boundsType)

    # Plots evolution of bridge distributions

    # Inputs:   - pInd is index of parameter
    #           - (N x np x Nϕ) array of particles
    #           - (N x Nϕ) matrix of weights
    #           - nLines is number of bridge distributions to plot
    #           - boundsType is either 0 (whole real line), 1 (positive real line), or 2 ((0,1) interval)



    Nϕ          = size(aParticles)[3]

    vn          = Int.( round.( collect(range(1,Nϕ,length = nLines)) ))

    transpStart = 0.5
    transpEnd   = 1
    vTransp     = collect(range(transpStart,transpEnd,length=nLines))

    k           = kde(aParticles[:,pInd,1], weights = Weights(mWeights[:,end]))
    st, en      = fFindPlotBounds(boundsType,k.x)


    vPlots      = Vector{Any}(undef,nLines)
    vPlots[1]   = plot(k.x[st:en],k.density[st:en],line=(:black,transpStart,2,:line),label="",xlabel=yLab,yticks=false,yaxis=false, xtickfont=font(16),ytickfont=font(16),guidefont=font(16),legendfontsize=16 )


    for nn in 2:nLines

        n       = vn[nn]

        k       = kde(aParticles[:,pInd,n], weights = Weights(mWeights[:,end]))
        st, en  = fFindPlotBounds(boundsType,k.x)

        vPlots[nn]   = plot(k.x[st:en],k.density[st:en],line=(:black,vTransp[nn],2,:line),label="",xlabel="",yticks=false,xticks=false,yaxis=false,xaxis=false)

    end

    yLims = ( minimum(getindex.(ylims.(vPlots),1)) , maximum(getindex.(ylims.(vPlots),2)) )
    xLims = ( minimum(getindex.(xlims.(vPlots),1)) , maximum(getindex.(xlims.(vPlots),2)) )

    return vPlots, xLims, yLims

end

function fPlotDistEvolGR(pInd,yLab,aParticles,mWeights,nLines,boundsType,xAxisTicks)

    # Plots evolution of bridge distributions

    # Inputs:   - pInd is index of parameter
    #           - (N x np x Nϕ) array of particles
    #           - (N x Nϕ) matrix of weights
    #           - nLines is number of bridge distributions to plot
    #           - boundsType is either 0 (whole real line), 1 (positive real line), or 2 ((0,1) interval)



    Nϕ          = size(aParticles)[3]

    vn          = Int.( round.( collect(range(1,Nϕ,length = nLines)) ))

    transpStart = 0.5
    transpEnd   = 1
    vTransp     = collect(range(transpStart,transpEnd,length=nLines))

    k           = kde(aParticles[:,pInd,1], weights = Weights(mWeights[:,end]))
    st, en      = fFindPlotBounds(boundsType,k.x)


    vx          = k.x[st:en] #range(-0.2,1.0,length=2048) #
    vy          = k.density[st:en] #range(0.0,2.0,length=2048) #
    ppp         = plot(ones(en-st+1)*nLines,vx,vy, camera=(75,45), line=(:black,transpStart,2,:line),label="",xlabel="n",ylabel=yLab,ylims=(k.x[st],k.x[en]), xtickfont=font(16),ytickfont=font(16),guidefont=font(16),legendfontsize=16 )
    # xflip!(true)


    for nn in 2:nLines

        n       = vn[nn]

        k       = kde(aParticles[:,pInd,n], weights = Weights(mWeights[:,end]))
        st, en  = fFindPlotBounds(boundsType,k.x)

        if nn < nLines
            ppp     = plot!(ones(en-st+1)*(Nϕ-n+1),k.x[st:en], k.density[st:en],camera=(75,45),line=(:black,vTransp[nn],2,:line),label="")
        else
            ppp     = plot!(ones(en-st+1)*(Nϕ-n+1),k.x[st:en], k.density[st:en],camera=(75,45),line=(:black,vTransp[nn],2,:line),label="",xticks=(vn,reverse(xAxisTicks)))
        end

    end

    return ppp

end


function fPlotDistEvol2D(pInd,yLab,aParticles,mWeights,nLines,boundsType,vLineLabs)

    # Plots evolution of bridge distributions

    # Inputs:   - pInd is index of parameter
    #           - (N x np x Nϕ) array of particles
    #           - (N x Nϕ) matrix of weights
    #           - nLines is number of bridge distributions to plot
    #           - boundsType is either 0 (whole real line), 1 (positive real line), or 2 ((0,1) interval)



    Nϕ          = size(aParticles)[3]

    vn          = Int.( round.( collect(range(1,Nϕ,length = nLines)) ))

    transpStart = 0.5
    transpEnd   = 1
    vTransp     = collect(range(transpStart,transpEnd,length=nLines))
    vThickness  = collect(range(1,2,length=nLines))

    vCols       = reverse(cgrad(:greys))

    k           = kde(aParticles[:,pInd,1], weights = Weights(mWeights[:,end]))
    st, en      = fFindPlotBounds(boundsType,k.x)


    vx          = k.x[st:en] #range(-0.2,1.0,length=2048) #
    vy          = k.density[st:en] #range(0.0,2.0,length=2048) #
    ppp         = plot(vx,vy,line=(vCols[0.4+0.5*1/nLines],transpStart,vThickness[1],:line), label=vLineLabs[1],xlabel=yLab,xlims=(k.x[st],k.x[en]),background_color_legend = nothing, xtickfont=font(16),ytickfont=font(16),guidefont=font(16),legendfontsize=16 )
    # xflip!(true)


    for nn in 2:nLines

        n       = vn[nn]

        k       = kde(aParticles[:,pInd,n], weights = Weights(mWeights[:,end]))
        st, en  = fFindPlotBounds(boundsType,k.x)

        ppp     = plot!(k.x[st:en], k.density[st:en],line=(vCols[0.4+0.5*nn/nLines],vTransp[nn],vThickness[nn],:line), label=vLineLabs[nn])

        # ppp     = plot!(ones(en-st+1)*n,k.x[st:en], k.density[st:en],camera=(60,30),line=(:black,vTransp[nn],2,:line),label="")

    end


    return ppp

end
