#!/usr/bin/env julia
include("Common.jl")
using Common
using RealNeuralNetworks.SWCs 
using HDF5
#using Plots
# use pyplot backend for vectorized output
#Plots.pyplot()

function main()
    args = parse_commandline()
    lengthList = Vector{Float64}()
    @sync begin
        for fileName in readdir( args["swcdir"] )
            if contains(fileName, ".swc")
                @async begin 
                    @time swc = SWCs.load( joinpath( args["swcdir"], fileName ) )
                    total_length = SWCs.get_total_length( swc )
                    println("length of $(fileName): $(total_length/1000) micron")
                    push!(lengthList, total_length)
                end
            end 
        end 
    end
    @show lengthList
    h5write("statistics.h5", "lengthList", lengthList)
    #p = histogram( lengthList )    
    # display(p)
end

@time main()
