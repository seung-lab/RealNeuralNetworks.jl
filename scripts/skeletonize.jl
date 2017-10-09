#!/usr/bin/env julia 
using ArgParse

using TEASAR
using TEASAR.NodeNets
using TEASAR.Manifests
using TEASAR.SWCs
using TEASAR.BranchNets

using JLD
using ProgressMeter

# this is the mip level 4
const EXPANSION = (16,16,1)
const MIP = 4

function trace(cellId; swcDir="/tmp/", jldDir="/tmp/", mip=MIP)
    manifest = Manifest("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4", 
                            "$(cellId):0", "gs://neuroglancer/zfish_v1/consensus-20170829/80_80_45")
    nodeNet = Manifests.trace(manifest, cellId)
    branchNet = BranchNet( nodeNet )
    swc = SWCs.SWC( BranchNet )
    SWCs.stretch_coordinates!(swc, mip)
    save(joinpath(jldDir, "$(cellId).jld"), "nodeNet", nodeNet, "branchNet", branchNet, "swc", swc)
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
function ArgParse.parse_item(::Type{UInt32}, x::AbstractString)
    return UInt32(parse(x))
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
        "--mip", "-m"
            help = "mip level of the dataset" 
            arg_type = UInt32
            default = UInt32(4)
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
        @showprogress 1 "tracing ..." for id in idList 
            trace(id; swcDir=args["swcdir"], jldDir=args["jlddir"], 
                        mip=args["mip"])
        end 
    else 
        trace(args["neuronid"]; swcDir = args["swcdir"], jldDir=args["jlddir"], 
                    mip=args["mip"])
    end 
end

main()
