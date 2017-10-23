#!/usr/bin/env julia
include("Common.jl")
using Common
using RealNeuralNetworks.SWCs 
using RealNeuralNetworks.BranchNets
using HDF5
#using Plots
# use pyplot backend for vectorized output
#Plots.pyplot()

function main()
    args = parse_commandline()
    numBranchingPointsList = Vector{Int}()
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

                    branchNet = BranchNet(swc)
                    numBranchingPoints = BranchNets.get_num_branching_points(branchNet)
                    println("number of branching points of cell $(fileName): $numBranchingPoints")
                    push!(numBranchingPointsList, numBranchingPoints)
                end
            end 
        end 
    end
    @show lengthList
    h5write("statistics.h5",    "totalPathLengthList",      totalPathLengthList,
                                "numBranchingPointsList",   numBranchingPointsList)
    #p = histogram( lengthList )    
    # display(p)
end

@time main()
