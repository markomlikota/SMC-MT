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
    xtickfont=font(14),ytickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14,background_color_legend = nothing
    )

    return p

end
