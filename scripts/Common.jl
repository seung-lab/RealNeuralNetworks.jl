module Common 

using ArgParse 
using JSON 

export parse_commandline 
export SEGMENT_ID, MIP, VOXEL_SIZE, SEGMENTATION_LAYER 
const MIP = UInt32(3)
const VOXEL_SIZE = (5,5,45)

function __init__()
    # setup AWS secrets 
    if isfile("/secrets/aws-secret.json")
        d = JSON.parsefile("/secrets/aws-secret.json")
        for (k,v) in d 
            ENV[k] = v 
        end 
    end 
end 

"""
customized argument parse for tuple
"""
function ArgParse.parse_item(::Type{NTuple{3,Int}}, x::AbstractString)
    return (map(Meta.parse, split(x,","))...,)
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
            arg_type = NTuple{3,Int}
            default = VOXEL_SIZE 
        "--mip", "-m"
            help = "mip level of the dataset" 
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

end # module
