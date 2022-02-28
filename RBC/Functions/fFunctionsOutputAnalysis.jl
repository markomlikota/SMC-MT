# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Dec 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Functions for output analysis



# -------------------------------------------------------------------------------


fGetSpecFolder(DGPspec,Modelspec,Priorspec,SMCspec) = string("dgp",DGPspec,"_m",Modelspec,"_pr",Priorspec,"_mc",SMCspec,"/")

fGetEstSpecFolderLT(ϕ_Nϕ) = string( "LTphiLast" , Int(ϕ_Nϕ * 100) , "/" )

fGetEstSpecFolderMT(ϕ_Nϕ_M0,Model0spec,Prior0spec,SMC0spec) = string( "MTphiLast" , Int(ϕ_Nϕ_M0 * 100) , "_m", Model0spec, "_pr", Prior0spec, "_mc", SMC0spec, "/" )

fGetEstSpecFolderMT(ϕ_Nϕ_M0,Model0spec,Prior0spec) = string( "MTphiLast" , Int(ϕ_Nϕ_M0 * 100) , "_m", Model0spec, "_pr", Prior0spec, "/" )


function fGetFileAllEstSpecs(DGPspec,SMCspec,nRun,sFileSuffix,vPhiLast)

    # Function to get some filetype for a given DGP and SMC specification and all estimation specifications (M1-LT, M1-MT, & corresponding M0-LT, where M1 and M0 defined below) and all runs.

    # Inputs:
    # - DGP and SMC specification files determining setup
    # - number of runs
    # - filename (e.g. _FinalStats.csv)

    # Output is vvvaStats, a vector (3-element; M0-LT, M1-MT, M1-LT) of vectors (first two have nPhiLast elements, M1-LT has one) of vectors (nRun elements) of arrays (actual content of the file; could be matrix, could be vector)


    #sSetupName      = fGetSetupName(VARspec,Priorspec,SMCspec)

    vvvaStats        = Vector{Any}[Vector{Any}(undef,1) for ss=1:3] # for the 3 vvSpecs


    # Obtain file across runs for M1-LT:

    Model1spec      = "L2"
    Prior1spec      = "L22"
    ϕ_Nϕ            = 1.0

    sSpecFolder     = fGetSpecFolder(DGPspec,Model1spec,Prior1spec,SMCspec)

    sEstSpecFolder  = fGetEstSpecFolderLT(ϕ_Nϕ)

    vaFiles           = Vector{Any}(undef,nRun)

    for rr = 1:nRun

        sPathHere     = sMyPath * "Output/" * sSpecFolder * "run" * string(rr) * "/" * sEstSpecFolder
        sFileName     = sFileSuffix         #sModelName * "_" * sFileSuffix
        dfFile        = fGetFile(sPathHere,sFileName)
        aFile         = convert(Matrix,dfFile)
        vaFiles[rr]   = aFile

    end

    vvvaStats[3][1] = vaFiles


    # Obtain file across runs for M1-MT and M0-LT:


    nPhiLast        = length(vPhiLast)

    for pphi in 1:nPhiLast

        ϕ_Nϕ_M0         = vPhiLast[pphi]
        #ϕ_Nϕ_M0 == 1.0 ? (vModelSpecIndices = 1:3) : (vModelSpecIndices = 1:2)


        # M0-LT:

        Model1spec      = "L1"
        Prior1spec      = "L12"
        ϕ_Nϕ            = deepcopy(ϕ_Nϕ_M0)

        sSpecFolder     = fGetSpecFolder(DGPspec,Model1spec,Prior1spec,SMCspec)

        sEstSpecFolder  = fGetEstSpecFolderLT(ϕ_Nϕ)

        vaFiles           = Vector{Any}(undef,nRun)

        for rr = 1:nRun

            sPathHere     = sMyPath * "Output/" * sSpecFolder * "run" * string(rr) * "/" * sEstSpecFolder
            sFileName     = sFileSuffix         #sModelName * "_" * sFileSuffix
            dfFile        = fGetFile(sPathHere,sFileName)
            aFile         = convert(Matrix,dfFile)
            vaFiles[rr]   = aFile

        end

        pphi == 1 ? vvvaStats[1][1] = vaFiles : push!(vvvaStats[1],vaFiles)


        # M1-MT:

        Model1spec      = "L2"
        Prior1spec      = "L22"
        Model0spec      = "L1"
        Prior0spec      = "L12"

        sSpecFolder     = fGetSpecFolder(DGPspec,Model1spec,Prior1spec,SMCspec)

        sEstSpecFolder  = fGetEstSpecFolderMT(ϕ_Nϕ_M0,Model0spec,Prior0spec)

        vaFiles           = Vector{Any}(undef,nRun)

        for rr = 1:nRun

            sPathHere     = sMyPath * "Output/" * sSpecFolder * "run" * string(rr) * "/" * sEstSpecFolder
            sFileName     = sFileSuffix         #sModelName * "_" * sFileSuffix
            dfFile        = fGetFile(sPathHere,sFileName)
            aFile         = convert(Matrix,dfFile)
            vaFiles[rr]   = aFile

        end

        pphi == 1 ? vvvaStats[2][1] = vaFiles : push!(vvvaStats[2],vaFiles)

    end


    return vvvaStats


end

function fGetStatIn_vvvaObject(vvvaObject,vPhiLast,ind,procedure,nosum=false)

    nPhilast = length(vPhiLast)

    if nosum == true
        vvObj = [getindex.(vvvaStats[3][1],ind),[getindex.(vvvaStats[2][pphi],ind) for pphi = 1:nPhiLast]...]
    else
        vvObj = [getindex.(vvvaStats[3][1],ind),[getindex.(vvvaStats[2][pphi],ind) + getindex.(vvvaStats[1][pphi],ind) for pphi = 1:nPhiLast]...]
    end

    if procedure == "mean"
        vStat   = mean.(vvObj)
    elseif procedure == "std"
        vStat   = std.(vvObj)
    elseif procedure == "lower"
        vStat   = map( vObj -> sort(vObj)[Int(nRun*0.05)], vvObj)
    elseif procedure == "upper"
        vStat   = map( vObj -> sort(vObj)[Int(nRun*0.95)], vvObj)
    end

    return vStat

end

function fMyRibbonPlotHere(vPhiLast,vy,vyDiff,sxLabel="",syLabel="")

    vxAxis      = [0.0,vPhiLast...]

    myBlue      = cgrad(:blues)[1.0]

    ppp         = plot(vxAxis,vy,line=(myBlue,0.9,3),label="",xlabel=sxLabel,title="",ylabel=syLabel,
    ribbon=vyDiff,fillalpha=0.2,fillcolor=myBlue, xtickfont=font(14),ytickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14)

    xticks!(vxAxis)

    return ppp

end

function fGetMedianTS(vaStageAllStats)

    # Computes vector of median (across runs) values for tempering parameter, one value for each stage of SMC.

    # Input: vector (nRun elements) of arrays with content from the file with suffix "_StageAll_Stats.csv"(output file of fSMC; in first column are values of tempering parameter across SMC-stages)

    nLongest    = maximum(size.(vaStageAllStats,1))
    nRun        = length(vaStageAllStats)
    vMedianTS   = zeros(nLongest)
    mTS         = Array{Float64,2}(undef,nLongest,nRun)
    fill!(mTS,NaN)

    for rr = 1:nRun

        nrow,ncol       = size(vaStageAllStats[rr])
        mTS[1:nrow,rr]  = vaStageAllStats[rr][:,1]

    end

    #vMedianTS   = [median(filter(!isnan,mTS[ii,:])) for ii=1:nLongest]
    vMedianTS   = mTS[:,1] # just take one realization (e.g. run number 1 here)

    return vMedianTS

end


function fPlotTSs(vPhiLastHere,vvvaStageAllStats,sCaption)

    # Function to create tempering-schedule plot.

    # Inputs:
    # - vector of phiLasts (i.e. \phi_{N_\phi}) to consider
    # - vector (3-element; M0-LT, M1-MT, M1-LT) of vectors (first two have nPhiLast elements, M1-LT has one) of vectors (nRun elements) of arrays with content from the file with suffix "StageAll_Stats.csv" (output file of fSMC; in first column are values of tempering parameter across SMC-stages)

    vCols       = cgrad(:blues)
    npphi       = length(vPhiLastHere)

    vTS_M1LT    = fGetMedianTS(vvvaStageAllStats[3][1])
    colHere     = vCols[1/(npphi+1)]
    TSsPlot     = plot(0:length(vTS_M1LT)-1,vTS_M1LT,line=(colHere,0.9,2,:dash),label="LT", legend=:bottomright,xlabel="n",title=sCaption,xtickfont=font(14),ytickfont=font(14),xguidefontsize=14,yguidefontsize=14,legendfontsize=14,background_color_legend = nothing)

    for pphi = 1:npphi

        phiLastHere     = vPhiLastHere[pphi]
        phiLastIndHere  = findfirst(vPhiLast .== phiLastHere)
        vTS_M1MT_Here   = fGetMedianTS(vvvaStageAllStats[2][pphi])
        sLab            = string("\$\\psi_*=",phiLastHere,"\$")
        colHere         = vCols[pphi/(npphi+1)]
        plot!(0:length(vTS_M1MT_Here)-1,vTS_M1MT_Here,line=(colHere,0.9,2),label=sLab)

        pphi == npphi ? (return TSsPlot) : nothing
    end

end


function fUpdateLimits(oldlimits,newlimits)

    # Takes the "union" of the two limits (intervals on ℜ)

    lower = min(oldlimits[1],newlimits[1])
    upper = max(oldlimits[2],newlimits[2])

    return (lower,upper)

end
