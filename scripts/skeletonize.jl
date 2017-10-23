#!/usr/bin/env julia 
include("Common.jl")
using .Common 

@everywhere using SQSChannels
@everywhere using GSDicts

@everywhere using RealNeuralNetworks
@everywhere using RealNeuralNetworks.NodeNets
@everywhere using RealNeuralNetworks.Manifests
@everywhere using RealNeuralNetworks.SWCs
@everywhere using RealNeuralNetworks.BranchNets

@everywhere function trace(cellId::Integer; swcDir      ::AbstractString = "/tmp/", 
                                mip         ::Integer        = MIP, 
                                voxelSize   ::Union{Tuple,Vector} = VOXEL_SIZE,
                                segmentationLayer::AbstractString = SEGMENTATION_LAYER)

    manifest = Manifest(joinpath(segmentationLayer, "mesh_mip_$(mip)"), "$(cellId):0", 
                        joinpath(segmentationLayer, 
                            "$(2^mip*voxelSize[1])_$(2^mip*voxelSize[2])_$(voxelSize[3])"))
    
    nodeNet = Manifests.trace(manifest, cellId)
    
    # transform to physical coordinate system
    NodeNets.stretch_coordinates!( nodeNet, mip )
    NodeNets.stretch_coordinates!( nodeNet, voxelSize)
 
   
   
    # reconnect the broken pieces and reset root to the soma center
    branchNet = BranchNet( nodeNet )
    swc = SWCs.SWC( branchNet )
    SWCs.save(swc, joinpath(swcDir, "$(cellId).swc"))
    
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip)"))
    d_str  = GSDict(joinpath(segmentationLayer, "swc"); valueType=String)
    d_bin["$cellId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$(cellId).swc"] = String(swc)
end 

"""
extract the cell id list from a text file which was transformed from google sheet shared in zfish_reconstruct channel
"""
function read_cell_id_list( fileName::String )
    idList = Vector{UInt32}()
    open(fileName) do f 
        for line in eachline(f)
            try 
                id = parse(split(line, ":")[1])
                push!(idList, UInt32(id))
            catch err
                if line != "\n"
                    rethrow()
                end 
            end 
        end 
    end
    return idList
end 

function main()
    args = parse_commandline()
    @show args
    if args["idlistfile"] != nothing
        idList = read_cell_id_list(args["idlistfile"])
        pmap(id -> trace(id; swcDir=args["swcdir"], mip=args["mip"]), 
                                                                                        idList)
    elseif args["sqsqueue"] != nothing 
        sqsChannel = SQSChannel( args["sqsqueue"] )
        while true 
            handle, body = fetch(sqsChannel) 
            d = JSON.parse( body )
            try 
                trace(d["id"]; swcDir=d["swcdir"], mip=d["mip"])
            catch err 
                if isa(err, KeyError)
                    warn("not tracable cell id: $(d["id"])") 
                else 
                    rethrow()
                end 
            end 
            delete!(sqsChannel, handle)
        end 
    else 
        trace(args["neuronid"]; swcDir = args["swcdir"],  
              mip=args["mip"], voxelSize=args["voxelsize"], 
              segmentationLayer=args["segmentationlayer"])
    end 
end

main()
