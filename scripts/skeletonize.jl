#!/usr/bin/env julia 
module ArgParses  

using ArgParse 
using JSON 

export parse_commandline 
export SEGMENT_ID, MIP, VOXEL_SIZE, SEGMENTATION_LAYER 
const MIP = UInt32(4)
const VOXEL_SIZE = map(Float32, (5,5,45))

function __init__()
    # setup AWS secrets 
    secret_file_path = expanduser("~/.cloudvolume/secrets/aws-secret.json")
    if !haskey(ENV, "AWS_ACCESS_KEY_ID") && isfile(secret_file_path)
        d = JSON.parsefile(secret_file_path)
        for (k,v) in d 
            ENV[k] = v 
        end 
    end 
end 

"""
customized argument parse for tuple
"""
function ArgParse.parse_item(::Type{NTuple{3,T}}, x::AbstractString) where T
    return (map(x->T(Meta.parse(x)), split(x,","))...,)
end
function ArgParse.parse_item(::Type{UInt32}, x::AbstractString)
    return UInt32(Meta.parse(x))
end 

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--neuronid", "-i"
            help = "the segment id to skeletonize"
            arg_type = Int
            # default = SEGMENT_ID #77497
        "--swcdir", "-s"
            help = "the directory to store swc file"
            arg_type = String
            default = "/tmp"
        "--outputpath", "-o"
            help = "the path to store outputs"
            arg_type = String
            default = "/tmp"
        "--voxelsize", "-v"
            help = "voxel size of the raw image, mip level 0"
            arg_type = NTuple{3,Float32}
            default = VOXEL_SIZE 
        "--mip", "-m"
            help = "mip level of the dataset. Note that our mip level start from 1." 
            arg_type = UInt32
            default = MIP 
        "--meshname", "-e"
            help = "folder name of mesh files"
            arg_type = String
            default = "mesh_mip_4"
        "--idlistfile", "-f"
            help = "the id list in text from google spreasheet"
            arg_type = String
        "--sqsqueue", "-q"
            help = "AWS SQS queue name"
            arg_type = String 
        "--segmentationlayer", "-l"
            help = "segmentation layer path in the cloud storage, only support Google Cloud Storage now"
            arg_type = String
    end 
    return parse_args(s)
end 

end # end of module ArgParses 

using .ArgParses 
using Distributed 

@everywhere using AWSSDK.SQS 
@everywhere using BigArrays.GSDicts

@everywhere using RealNeuralNetworks
@everywhere using RealNeuralNetworks.NodeNets
@everywhere using RealNeuralNetworks.Manifests
@everywhere using RealNeuralNetworks.SWCs
@everywhere using RealNeuralNetworks.Neurons


@everywhere function trace(neuronId::Integer; swcDir      ::AbstractString = "/tmp/", 
                                mip         ::Integer        = MIP, 
                                meshName    ::String = "mesh_mip_$MIP",
                                voxelSize   ::Union{Tuple,Vector} = VOXEL_SIZE,
                                segmentationLayer::AbstractString = SEGMENTATION_LAYER)
    println("skeletonize neuron: $(neuronId)")
    println("fetching manifest...")
    manifest = Manifest(joinpath(segmentationLayer, meshName), "$(neuronId):0", 
                        joinpath(segmentationLayer), mip)
    
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
    neuron = Neurons.remove_subtree_in_soma(neuron)
    neuron = Neurons.remove_hair( neuron )

    swc = SWCs.SWC( neuron )
    SWCs.save(swc, joinpath(swcDir, "$(neuronId).swc"))
    
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip-1)"))
    d_str  = GSDict(joinpath(segmentationLayer, "swc"))
    d_swc_bin  = GSDict(joinpath(segmentationLayer, "swc_bin"))
    d_bin["$neuronId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$(neuronId).swc"] = String(swc)
    d_swc_bin["$(neuronId).swc.bin"] = SWCs.serialize(swc)
    
    println("postprocessing the neuron...")
    neuron = Neurons.postprocessing(neuron)
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip-1)_postprocessed"))
    d_str  = GSDict(joinpath(segmentationLayer, "postprocessed_swc"))
    d_swc_bin  = GSDict(joinpath(segmentationLayer, "postprocessed_swc_bin"))
    d_bin["$neuronId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$(neuronId).swc"] = String(swc)
    d_swc_bin["$(neuronId).swc.bin"] = SWCs.serialize(swc)

    println("resample the neuron...")
    neuron = Neurons.resample(neuron, Float32(500))
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip-1)_resampled"))
    d_str  = GSDict(joinpath(segmentationLayer, "resampled_swc"))
    d_swc_bin  = GSDict(joinpath(segmentationLayer, "resampled_swc_bin"))
    d_bin["$neuronId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$(neuronId).swc"] = String(swc)
    d_swc_bin["$(neuronId).swc.bin"] = SWCs.serialize(swc)

    println("finish skeletonize.")
end 

"""
extract the cell id list from a text file which was transformed from google sheet shared in zfish_reconstruct channel
"""
function read_cell_id_list( fileName::String )
    idList = Vector{UInt32}()
    open(fileName) do f 
        for line in eachline(f)
            try 
                id_list = map(Meta.parse, split(line, ","))
                for id in id_list 
                    push!(idList, UInt32(id))
                end
            catch err
                if line != "\n"
                    rethrow()
                end 
            end 
        end 
    end
    return idList
end

@everywhere function trace_queue(args::Dict{String, Any})
    queueName = args["sqsqueue"]
    queueUrl = SQS.get_queue_url(QueueName=queueName)["QueueUrl"]
    while true
        println("try to pull task from SQS queue: $(args["sqsqueue"])")
        local message 
        try 
            message = SQS.receive_message(QueueUrl=queueUrl)["messages"][1]
        catch err 
            @warn("no task fetched, could be all done!")
            println(err)
            break
            #println("sleep for 30 secs and then retry")
            #sleep(30)
            #continue 
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
        @sync for pid in 1:Distributed.nworkers()
            remote_do(trace_queue(args), pid)
        end
    else 
        trace(args["neuronid"]; swcDir = args["swcdir"],  
              mip=args["mip"], voxelSize=args["voxelsize"],
              meshName=args["meshname"],
              segmentationLayer=args["segmentationlayer"])
    end 
end

main()
