#!/usr/bin/env julia 
using ArgParse

using TEASAR
using TEASAR.Skeletons
using TEASAR.Manifests
using TEASAR.SWCs
using JLD
using ProgressMeter

const VOXEL_SIZE = (80,80,45)

function skeletonize(cellId; swcDir="/tmp/", jldDir="/tmp/", voxel_size=VOXEL_SIZE)
    manifest = Manifest("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4", 
                            "$(cellId):0", "gs://neuroglancer/zfish_v1/consensus-20170829/80_80_45")
    skeleton = Manifests.trace(manifest, cellId)
    swc = SWC( skeleton )
    save(joinpath(jldDir, "$(cellId).jld"), "skeleton", skeleton, "swc", swc)
    SWCs.stretch_coordinates!(swc, voxel_size)
    SWCs.save(swc, joinpath(swcDir, "$(cellId).swc"))
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

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--neuronid", "-i"
            help = "the segment id to skeletonize"
            arg_type = Int
            default = 77497
        "--swcdir", "-s"
            help = "the directory to store swc file"
            arg_type = String
            default = "/tmp"
        "--jlddir", "-j"
            help = "the directory to store jld file"
            arg_type = String
            default = "/tmp"
        "--voxelsize", "-v"
            help = "voxel size of dataset" 
            arg_type = NTuple{3, Int}
            default = VOXEL_SIZE
        "--idlistfile", "-f"
            help = "the id list in text from google spreasheet"
            arg_type = String
    end 
    return parse_args(s)
end 

function main()
    args = parse_commandline()
    @show args
    if args["idlistfile"] != nothing
        idList = read_cell_id_list(args["idlistfile"])
        @showprogress 1 "skeletonize ..." for id in idList 
            skeletonize(id; swcDir=args["swcdir"], jldDir=args["jlddir"], 
                        voxel_size=args["voxelsize"])
        end 
    else 
        skeletonize(args["neuronid"]; swcDir = args["swcdir"], jldDir=args["jlddir"], 
                    voxel_size=args["voxelsize"])
    end 
end

main()
