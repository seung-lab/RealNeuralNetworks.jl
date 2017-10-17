#!/usr/bin/env julia 
using ArgParse
using GSDicts

using RealNeuralNetworks
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.Manifests
using RealNeuralNetworks.SWCs
using RealNeuralNetworks.BranchNets

# this is the mip level 4
const SEGMENTAT_ID = 92540687
const MIP = UInt32(4)
const VOXEL_SIZE = (5,5,45)
const SEGMENTATION_LAYER ="gs://neuroglancer/zfish_v1/consensus-20170928"

function trace(cellId::Integer; swcDir      ::AbstractString = "/tmp/", 
                                jldDir      ::AbstractString = "/tmp/", 
                                mip         ::Integer        = MIP, 
                                voxelSize   ::Union{Tuple,Vector} = VOXEL_SIZE,
                                segmentationLayer::AbstractString = SEGMENTATION_LAYER)

    manifest = Manifest(joinpath(segmentationLayer, "mesh_mip_$(mip)"), "$(cellId):0", 
                        joinpath(segmentationLayer, 
                            "$(2^mip*voxelSize[1])_$(2^mip*voxelSize[2])_$(voxelSize[3])"))
    
    nodeNet = Manifests.trace(manifest, cellId)
    
    # transform to physical coordinate system
    NodeNets.add_offset!( nodeNet, (-1,-1,-1) )
    NodeNets.stretch_coordinates!( nodeNet, mip )
    NodeNets.stretch_coordinates!( nodeNet, voxelSize)
 
   
   
    # reconnect the broken pieces and reset root to the soma center
    branchNet = BranchNet( nodeNet )
    swc = SWCs.SWC( branchNet )
    SWCs.save(swc, joinpath(swcDir, "$(cellId).swc"))
    
    # save to neuroglancer
    d_bin  = GSDict(joinpath(segmentationLayer, "skeleton_mip_$(mip)"))
    d_str  = GSDict(joinpath(segmentationLayer, "swc"; valueType=String))
    d_bin["$cellId"] = SWCs.get_neuroglancer_precomputed( swc )
    d_str["$cellId"] = String(swc)
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

"""
customized argument parse for tuple
"""
function ArgParse.parse_item(::Type{NTuple{3,Int}}, x::AbstractString)
    return map(parse, split(x,","))
end
function ArgParse.parse_item(::Type{UInt32}, x::AbstractString)
    return UInt32(parse(x))
end 

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--neuronid", "-i"
            help = "the segment id to skeletonize"
            arg_type = Int
            default = SEGMENTAT_ID #77497
        "--swcdir", "-s"
            help = "the directory to store swc file"
            arg_type = String
            default = "/tmp"
        "--jlddir", "-j"
            help = "the directory to store jld file"
            arg_type = String
            default = "/tmp"
        "--voxelsize", "-v"
            help = "voxel size of the raw image, mip level 0"
            arg_type = NTuple{3,Int}
            default = VOXEL_SIZE 
        "--mip", "-m"
            help = "mip level of the dataset" 
            arg_type = UInt32
            default = MIP 
        "--idlistfile", "-f"
            help = "the id list in text from google spreasheet"
            arg_type = String
        "--segmentationlayer", "-l"
            help = "segmentation layer path in the cloud storage, only support Google Cloud Storage now"
            arg_type = String
            default = SEGMENTATION_LAYER
    end 
    return parse_args(s)
end 

function main()
    args = parse_commandline()
    @show args
    if args["idlistfile"] != nothing
        idList = read_cell_id_list(args["idlistfile"])
        for id in idList 
            trace(id; swcDir=args["swcdir"], jldDir=args["jlddir"], 
                        mip=args["mip"])
        end 
    else 
        trace(args["neuronid"]; swcDir = args["swcdir"], jldDir=args["jlddir"], 
              mip=args["mip"], voxelSize=args["voxelsize"], 
              segmentationLayer=args["segmentationlayer"])
    end 
end

main()
