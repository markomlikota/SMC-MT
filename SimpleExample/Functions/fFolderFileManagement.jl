# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# November 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Various functions to deal with files and folders for output storage.



# -------------------------------------------------------------------------------



function fAppendFolder(sPath,sFolderName)

    # Searches for last occurence of a folder starting with "sFolderName" (and followed by some integer) in "sPath", and creates a folder starting with "sFolderName" and the highest integer + 1. Returns this sFolderName.
    # e.g. "SetOfRuns1" and "SetOfRuns2" exist in "sPath", function creates and returns "SetOfRuns3".

    vFileNames = readdir(sPath)

    vInd       = findall(occursin.(sFolderName,vFileNames) .== 1)

    if length(vInd) == 0
        folder   = string(sFolderName,1)
    else
        maxInd = maximum( [ parse(Int,replace(string(vFileNames[vInd][i]),sFolderName=>"")) for i = 1:length(vInd) ] )
        folder   = string(sFolderName,maxInd+1)
    end

    try
        mkdir(string(sPath,"/",folder))
    catch
    end

    return folder

end


function fFindNewestFolder(sPath,sFolderName)

    # Searches for last (newest) occurence of a folder starting with "sFolderName" (and followed by some integer) in "sPath", and returns this foldername and its number (e.g. "SetOfRuns3", 3)

    vFolderNames = readdir(sPath)

    vInd         = findall(occursin.(sFolderName,vFolderNames) .== 1)

    if length(vInd) == 0
        return nothing
    else
        maxInd = maximum( [ parse(Int,replace(string(vFolderNames[vInd][i]),sFolderName=>"")) for i = 1:length(vInd) ] )
        folder   = string(sFolderName,maxInd)
        return folder, maxInd
    end

end


function fFindNewestFile(sPath,sFileName,sFileType)

    # Searches for last (newest) occurence of a file starting with "sFileName" (and followed by some integer) in "sPath", and returns this file name and its number (e.g. "PostDraws4.csv", 3)
    # sFileType is ".csv" for example

    vFileNames = readdir(sPath)

    vInd       = findall(occursin.(sFileName,vFileNames) .== 1)

    if length(vInd) == 0
        return nothing
    else
        maxInd = maximum( [ parse(Int,split(replace(string(vFileNames[vInd][i]),sFileName=>""),".")[1]) for i = 1:length(vInd) ] )
        file   = string(sFileName,maxInd,sFileType)
        return file, maxInd
    end

end


function fCombinevStrings(vStrings,sInBetw="_")

    # Input: vStrings is vector of strings
    # Output: string composed of combined elements in args, separated by underlines in-between

    vStrings_       = string.(vStrings,sInBetw)

    vStrings_[end]  = split(vStrings_[end],sInBetw)[1]

    return string(vStrings_...)

end


function fGetFile(sPath,sFileName)

    f                   = open(sPath * sFileName)
    mC                  = readdlm(f, ',',skipstart=1) #read in CSV file

    return mC

end


function fResetCSV(sFileName)

    # Deletes all contents in csv-file "filename" (writes empty data frame to it)

    CSV.write( sFileName, DataFrame() )

end


function fPrintDF(df,nDigits)

    # Prints data frame df as string with specified number of digits

    nR, nC                      = size(df)

    vIsFloat                    = eltype.(eachcol(df))    .== Float64

    vIndFloatCols               = collect(1:nC)[vIsFloat]

    [df[!,propertynames(df)[cc]] = map( (x) -> string(round(x; digits=nDigits)) , df[!,propertynames(df)[cc]]) for cc in vIndFloatCols]

    for cc in vIndFloatCols, rr in 1:nR

        nDigs       = length(split(df[rr,cc],".")[2])
        nDigs < nDigits ? df[rr,cc] = string(df[rr,cc],string(repeat(["0"],nDigits-nDigs)...)) : nothing

    end

    return df

end
