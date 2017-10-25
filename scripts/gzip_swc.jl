#!/usr/bin/env julia
using Libz
using RealNeuralNetworks
using RealNeuralNetworks.SWCs
include("Common.jl"); using .Common 

function main()
    args = parse_commandline()
    for fileName in readdir( args["swcdir"] )
        println("transfering file: $(fileName)")
        swcString = readstring( joinpath(args["swcdir"], fileName) )
        outputFileName = joinpath(args["outputpath"], split(fileName, ".")[1] * ".swc.gz")
        io = open(outputFileName, "w")
        stream = ZlibDeflateOutputStream(io)
        write(stream, swcString)
        close(stream)
    end 
end 

main()
