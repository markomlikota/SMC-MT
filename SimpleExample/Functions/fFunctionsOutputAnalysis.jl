# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Various functions for output analysis codes.



# -------------------------------------------------------------------------------



function fGetPlotOneDistanceMeasure(mParamInfo,mDistanceMeasures,vStdIndicesToUse,measureInd)

    fvIndicesForStd(x) = mParamInfo[:,4] .== x
    fLabelForStd(x) = string("σ0 = ",x)
    vMeasureLabsShort = ["Area under product of pdfs","Kullback-Leibler Distance","Area under minimum of pdfs","V-M_0 of pdf ratio"]
    vMeasureLabs = string.(vMeasureLabsShort," (target=N(0,1))")

    ss = 1

    ssInd         = vStdIndicesToUse[ss]

    myBlue        = cgrad(:blues)[ssInd/nStds0]

    σ             = vStds0[ssInd]

    pp = scatter(mParamInfo[fvIndicesForStd(σ),3],mDistanceMeasures[fvIndicesForStd(σ),measureInd],color=myBlue,markerstrokewidth=0.5,markersize=10,label=fLabelForStd(σ),xlabel="μ0",title=vMeasureLabs[measureInd])

    for ss = 2:length(vStdIndicesToUse)
      ssInd         = vStdIndicesToUse[ss]

      myBlue        = cgrad(:blues)[ssInd/nStds0]

      σ             = vStds0[ssInd]

      scatter!(mParamInfo[fvIndicesForStd(σ),3],mDistanceMeasures[fvIndicesForStd(σ),measureInd],color=myBlue,markerstrokewidth=0.5,markersize=10,label=fLabelForStd(σ))
    end

    return pp

end


function fGetPlotTwoDistanceMeasures(mParamInfo,mDistanceMeasures,vStdIndicesToUse,vMeasureInd)

    fvIndicesForStd(x) = mParamInfo[:,4] .== x
    fLabelForStd(x) = string("σ0 = ",x)
    vMeasureLabsShort = ["Area under product of pdfs","Kullback-Leibler Distance","Area under minimum of pdfs","V-M_0 of pdf ratio"]
    vMeasureLabs = string.(vMeasureLabsShort," (target=N(0,1))")

    ss = 1

    ssInd         = vStdIndicesToUse[ss]

    myBlue        = cgrad(:blues)[ssInd/nStds0]

    σ             = vStds0[ssInd]

    pp = scatter(mDistanceMeasures[fvIndicesForStd(σ),vMeasureInd[1]],mDistanceMeasures[fvIndicesForStd(σ),vMeasureInd[2]],color=myBlue,markerstrokewidth=0.5,markersize=10,label=fLabelForStd(σ),xlabel=vMeasureLabs[vMeasureInd[1]],ylabel=vMeasureLabs[vMeasureInd[2]])

    for ss = 2:length(vStdIndicesToUse)
      ssInd         = vStdIndicesToUse[ss]

      myBlue        = cgrad(:blues)[ssInd/nStds0]

      σ             = vStds0[ssInd]

      scatter!(mDistanceMeasures[fvIndicesForStd(σ),vMeasureInd[1]],mDistanceMeasures[fvIndicesForStd(σ),vMeasureInd[2]],color=myBlue,markerstrokewidth=0.5,markersize=10,label=fLabelForStd(σ))
    end

    return pp

end
