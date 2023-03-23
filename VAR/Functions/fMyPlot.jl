# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Simple function for plotting one vector (vy) against another (vx)



# -------------------------------------------------------------------------------


function fMyPlot(vx,vy,sCap="",sxLabel="",syLabel="",sLineLabel="",col=:black)

    # Standard plot function for line plot

    p       =
    plot(
    vx, vy,
    title=sCap, xlabel=sxLabel, ylabel=syLabel, label=sLineLabel,
    grid=false,line=(col,0.9,2,:line),
    xtickfont=font(16),ytickfont=font(16),guidefont=font(16),legendfontsize=16,background_color_legend = nothing
    )

    return p

end
