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
	sqsChannel = SQSChannel( args["sqsqueue"] )
	for id in idList
        @show id
        args["id"] = id
		put!(sqsChannel, JSON.json(args))
	end 
end

main()
