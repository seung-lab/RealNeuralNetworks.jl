#!/usr/bin/env julia 
include("Common.jl")
using .Common 

using Distributed 

@everywhere using AWSSDK.SQS 
@everywhere using BigArrays.GSDicts

@everywhere using RealNeuralNetworks
@everywhere using RealNeuralNetworks.NodeNets
@everywhere using RealNeuralNetworks.Manifests
@everywhere using RealNeuralNetworks.SWCs
@everywhere using RealNeuralNetworks.Neurons


if !haskey(ENV, "AWS_ACCESS_KEY_ID") && isfile("/secrets/aws-secret.json")
    d = JSON.parsefile("/secrets/aws-secret.json")
    for (k,v) in d
        ENV[k] = v
    end 
end 

@everywhere function trace(neuronId::Integer; swcDir      ::AbstractString = "/tmp/", 
                                mip         ::Integer        = MIP, 
                                meshName    ::String = "mesh_mip_$MIP",
                                voxelSize   ::Union{Tuple,Vector} = VOXEL_SIZE,
                                segmentationLayer::AbstractString = SEGMENTATION_LAYER)
    println("fetching manifest...")
    manifest = Manifest(joinpath(segmentationLayer, meshName), "$(neuronId):0", 
                        joinpath(segmentationLayer, 
                            "$(2^mip*voxelSize[1])_$(2^mip*voxelSize[2])_$(voxelSize[3])"))
    
    println("skeletonizing the point cloud...")
    nodeNet = Manifests.trace(manifest)
    # transform to physical coordinate system
    # the mesh coordinate start from 0 rather then 1
    NodeNets.add_offset!(nodeNet, (-one(Float32), -one(Float32), -one(Float32)))
    NodeNets.stretch_coordinates!( nodeNet, mip )
    NodeNets.stretch_coordinates!( nodeNet, voxelSize)
    NodeNets.save(nodeNet, "/tmp/$(neuronId).swc")
 
    # reconnect the broken pieces and reset root to the soma center 
    neuron = Neuron( nodeNet )
    #neuron = Neurons.remove_subtree_in_soma(neuron)
    #neuron = Neurons.remove_hair( neuron )

    swc = SWCs.SWC( neuron )
    SWCs.save(swc, joinpath(swcDir, "$(neuronId).swc"))
    
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip)"))
    d_str  = GSDict(joinpath(segmentationLayer, "swc"); valueType=String)
    d_swc_bin  = GSDict(joinpath(segmentationLayer, "swc_bin"))
    d_bin["$neuronId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$(neuronId).swc"] = String(swc)
    d_swc_bin["$(neuronId).swc.bin"] = SWCs.serialize(swc)
    
    println("postprocessing the neuron...")
    neuron = Neurons.postprocessing(neuron)
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip)_postprocessed"))
    d_str  = GSDict(joinpath(segmentationLayer, "postprocessed_swc"); valueType=String)
    d_swc_bin  = GSDict(joinpath(segmentationLayer, "postprocessed_swc_bin"))
    d_bin["$neuronId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$(neuronId).swc"] = String(swc)
    d_swc_bin["$(neuronId).swc.bin"] = SWCs.serialize(swc)

    println("resample the neuron...")
    neuron = Neurons.resample(neuron, Float32(500))
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip)_resampled"))
    d_str  = GSDict(joinpath(segmentationLayer, "resampled_swc"); valueType=String)
    d_swc_bin  = GSDict(joinpath(segmentationLayer, "resampled_swc_bin"))
    d_bin["$neuronId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$(neuronId).swc"] = String(swc)
    d_swc_bin["$(neuronId).swc.bin"] = SWCs.serialize(swc)

end 

"""
extract the cell id list from a text file which was transformed from google sheet shared in zfish_reconstruct channel
"""
function read_cell_id_list( fileName::String )
    idList = Vector{UInt32}()
    open(fileName) do f 
        for line in eachline(f)
            try 
                id = Meta.parse(split(line, ":")[1])
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
        pmap(id -> trace(id; swcDir=args["swcdir"], mip=args["mip"], 
                         meshName=args["meshname"],
                         voxelSize = args["voxelsize"],
                         segmentationLayer=args["segmentationlayer"]), idList)
    elseif args["sqsqueue"] != nothing 
        queueUrl = SQS.get_queue_url(QueueName=args["sqsqueue"])["QueueUrl"]
        while true
            println("try to pull task from SQS queue: $(args["sqsqueue"])")
            local message 
            try 
                message = SQS.receive_message(QueueUrl=queueUrl)["messages"][1]
            catch err 
                @warn("no task fetched, could be all done!")
                println(err)
                println("sleep for 30 secs and then retry")
                sleep(30)
                continue 
            end 
            receiptHandle = message["ReceiptHandle"]
            id = Meta.parse(message["Body"])
            println("tracing cell: $(id)")
            try 
                trace(id; swcDir=args["swcdir"], mip=args["mip"], 
                      meshName=args["meshname"],
                      voxelSize = args["voxelsize"],
                      segmentationLayer = args["segmentationlayer"])
            catch err 
                if isa(err, KeyError)
                    println("key not found: $(id)")
                else 
                    rethrow() 
                end 
            end
            println("delete task in queue: $(receiptHandle)")
            SQS.delete_message(QueueUrl=queueUrl, ReceiptHandle=receiptHandle)
        end 
    else 
        trace(args["neuronid"]; swcDir = args["swcdir"],  
              mip=args["mip"], voxelSize=args["voxelsize"],
              meshName=args["meshname"],
              segmentationLayer=args["segmentationlayer"])
    end 
end

main()
