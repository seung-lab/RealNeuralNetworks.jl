module Common

using ArgParse 
export parse_commandline 

# this is the mip level 4
const SEGMENTAT_ID = 92540687
const MIP = UInt32(4)
const VOXEL_SIZE = (5,5,45)
const SEGMENTATION_LAYER ="gs://neuroglancer/zfish_v1/consensus-20170928"

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
        "--segmentationlayer", "-l"
            help = "segmentation layer path in the cloud storage, only support Google Cloud Storage now"
            arg_type = String
            default = SEGMENTATION_LAYER
    end 
    return parse_args(s)
end 

end # module
