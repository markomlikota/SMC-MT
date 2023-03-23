# Mlikota & Schorfheide : Sequential Monte Carlo With Model Tempering
# Nov 2021, mlikota@sas.upenn.edu

# -------------------------------------------------------------------------------


# Various functions to deal with files and folders for output storage.


# Note: replace DataFrame(x) with DataFrame(x,:auto) in CSV.write once newest version of CSV package works with Julia 1.6.3 !
# Note: replace readdlm(x, ',', Float64, '\n'; skipstart=1) with CSV.read(x;header=1) once newest version of CSV package works with Julia 1.6.3 !
# i.e. change definitions of fMyCSVREAD, fMyCSVWRITE and fMyCSVAPPEND on bottom of this file


# -------------------------------------------------------------------------------

# For Marko's M1-Macbook: (newest version of CSV package doesn't work on new Mac-chips apparently)


fMyCSVREAD(sPathToRead) = readdlm(sPathToRead, ',', Float64, '\n'; skipstart=1)

function fMyCSVWRITE(sPathToWrite,mx,vHeader=0)

    if vHeader != 0
        CSV.write( sPathToWrite,  DataFrame(permutedims(vHeader)), header=false )
        CSV.write( sPathToWrite,  DataFrame(mx), append=true )
    else
        CSV.write( sPathToWrite, DataFrame(mx) )
    end

end

fMyCSVAPPEND(sPathToWrite,mx) = CSV.write( sPathToWrite,  DataFrame(mx), append=true )



# -------------------------------------------------------------------------------

# For Frank:


# fMyCSVREAD(sPathToRead) = CSV.read(sPathToRead;header=1)
#
# function fMyCSVWRITE(sPathToWrite,mx,vHeader=0)
#
#     if vHeader != 0
#         CSV.write( sPathToWrite,  DataFrame(permutedims(vHeader),:auto), header=false )
#         CSV.write( sPathToWrite,  DataFrame(mx,:auto), append=true )
#     else
#         CSV.write( sPathToWrite, DataFrame(mx,:auto) )
#     end
#
# end
#
# fMyCSVAPPEND(sPathToWrite,mx) = CSV.write( sPathToWrite,  DataFrame(mx,:auto), append=true )
