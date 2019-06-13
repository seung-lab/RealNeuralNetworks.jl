#!/usr/bin/env julia
include("Common.jl")
using Common
using RealNeuralNetworks.SWCs 
using RealNeuralNetworks.Neurons
using HDF5
#using Plots
# use pyplot backend for vectorized output
#Plots.pyplot()

function main()
    args = parse_commandline()
    numSegmentingPointsList = Vector{Int}()
    totalPathLengthList = Vector{Float64}()
    @sync begin
        for fileName in readdir( args["swcdir"] )
            if contains(fileName, ".swc")
                @async begin 
                    @time swc = SWCs.load( joinpath( args["swcdir"], fileName ) )
                    totalPathLength = SWCs.get_total_length( swc )
                    println("total path length of cell $(fileName): " * 
                            "$(totalPathLength/1000) micron")
                    push!(totalPathLengthList, totalPathLength)

                    neuron = Neuron(swc)
                    numSegmentingPoints = Neurons.get_num_branching_points(neuron)
                    println("number of branching points of cell $(fileName): $numSegmentingPoints")
                    push!(numSegmentingPointsList, numSegmentingPoints)
                end
            end 
        end 
    end
    h5write("statistics.h5",    "totalPathLengthList",      totalPathLengthList)
    h5write("statistics.h5",    "numSegmentingPointsList",   numSegmentingPointsList)
    #p = histogram( lengthList )    
    # display(p)
end

@time main()
