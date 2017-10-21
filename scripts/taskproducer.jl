#!/usr/bin/env julia
using SQSChannels 
using ArgParse 
using JSON

# this is the mip level 4
@everywhere const SEGMENT_ID = 92540687
@everywhere const MIP = UInt32(4)
@everywhere const VOXEL_SIZE = (5,5,45)
@everywhere const SEGMENTATION_LAYER ="gs://neuroglancer/zfish_v1/consensus-20170928"


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--neuronid", "-i"
            help = "the segment id to skeletonize"
            arg_type = Int
            default = SEGMENT_ID #77497
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
        "--sqsqueue", "-q"
            help = "AWS SQS queue name"
            arg_type = String
            default = "skeleton"
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
    @assert args["idlistfile"] != nothing
     
    idList = map( parse, split( readstring(args["idlistfile"]), "\n" ) )
	sqsChannel = SQSChannel( args["sqsqueue"] )
	for id in idList
        @show id
        args["id"] = id
		put!(sqsChannel, JSON.json(args))
	end 
end

main()
