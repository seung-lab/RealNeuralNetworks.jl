#!/usr/bin/env julia 
using ArgParse

using TEASAR
using TEASAR.Skeletons
using TEASAR.Manifests
using TEASAR.SWCs
using JLD

const VOXEL_SIZE = (80,80,45)

function skeletonize(cellId; swcDir="/tmp/", voxel_size=VOXEL_SIZE)
    manifest = Manifest("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4", 
                            "$(cellId):0", "gs://neuroglancer/zfish_v1/consensus-20170829/80_80_45")
    skeleton = Manifests.trace(manifest, cellId)
    swc = SWC( skeleton )
    save("/tmp/$(cellId).jld", "skeleton", skeleton, "swc", swc)
    SWCs.stretch_coordinates!(swc, voxel_size)
    SWCs.save(swc, joinpath(swcDir, "$(cellId).swc"))
end 

"""
customized argument parse for tuple
"""
function ArgParse.parse_item(::Type{NTuple{3,Int}}, x::AbstractString)
    return map(parse, split(x,","))
end 

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--neuronid", "-i"
            help = "the segment id to skeletonize"
            arg_type = Int
            default = 77497
        "--swcdir", "-d"
            help = "the directory to store swc file"
            arg_type = String
            default = "/tmp"
        "--voxelsize", "-v"
            help = "voxel size of dataset" 
            arg_type = NTuple{3, Int}
            default = VOXEL_SIZE 
    end 
    return parse_args(s)
end 

function main()
    args = parse_commandline()
    @show args 
    skeletonize(args["neuronid"]; swcDir = args["swcdir"], voxel_size=args["voxelsize"])
end

main()
