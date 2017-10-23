#!/usr/bin/env julia
include("Common.jl")
using .Common 
using SQSChannels 
using JSON

function main()
    args = parse_commandline()
    @show args
    @assert args["idlistfile"] != nothing
     
    idList = map( parse, split( readstring(args["idlistfile"]), "\n" ) )
    argList = Vector{String}()
    for id in idList 
        args["id"] = id 
        push!(argList, JSON.json(args))
    end

    # put to AWS SQS queue
	sqsChannel = SQSChannel( args["sqsqueue"] )
    put!(sqsChannel, argList )
end

main()
