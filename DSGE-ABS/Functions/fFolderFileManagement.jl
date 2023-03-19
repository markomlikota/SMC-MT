# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Various functions to deal with files and folders for output storage.


# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !
# i.e. change definitions of fMyCSVREAD, fMyCSVWRITE and fMyCSVAPPEND on bottom of this file


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

function fDoesFileExist(sPath,sFileName)

    vFileNames = readdir(sPath)

    vInd       = findall(occursin.(sFileName,vFileNames) .== 1)

    return length(vInd) != 0

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



# -------------------------------------------------------------------------------

# For Marko's M1-Macbook: (newest version of CSV package doesn't work on new Mac-chips apparently)


# fMyCSVREAD(sPathToRead) = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)
#
# function fMyCSVWRITE(sPathToWrite,mx,vHeader=0)
#
#     if vHeader != 0
#         CSV.write( sPathToWrite,  DataFrame(permutedims(vHeader)), header=false )
#         CSV.write( sPathToWrite,  DataFrame(mx), append=true )
#     else
#         CSV.write( sPathToWrite, DataFrame(mx) )
#     end
#
# end
#
# fMyCSVAPPEND(sPathToWrite,mx) = CSV.write( sPathToWrite,  DataFrame(mx), append=true )
#


# -------------------------------------------------------------------------------

# For Frank:


fMyCSVREAD(sPathToRead) = CSV.read(sPathToRead;header=1)

function fMyCSVWRITE(sPathToWrite,mx,vHeader=0)

    if vHeader != 0
        CSV.write( sPathToWrite,  DataFrame(permutedims(vHeader),:auto), header=false )
        CSV.write( sPathToWrite,  DataFrame(mx,:auto), append=true )
    else
        CSV.write( sPathToWrite, DataFrame(mx,:auto) )
    end

end

fMyCSVAPPEND(sPathToWrite,mx) = CSV.write( sPathToWrite,  DataFrame(mx,:auto), append=true )
