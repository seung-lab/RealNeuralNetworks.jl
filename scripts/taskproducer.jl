#!/usr/bin/env julia
include("Common.jl")
using .Common 
using JSON
using AWSSDK.SQS 
using AWSCore
using ProgressMeter

const AWS_CREDENTIAL = AWSCore.aws_config()

function main()
    args = parse_commandline()
    @show args
    @assert args["idlistfile"] != nothing
     
    idList = map( parse, split( readstring(args["idlistfile"]), "\n" ) )
    messageList = Vector{String}()
    for id in idList 
        push!(messageList, string(id))
    end

    # put to AWS SQS queue
    const queueUrl = SQS.get_queue_url( QueueName=args["sqsqueue"] )["QueueUrl"]

#    for message in messageList
#        @show message
#        SQS.send_message( QueueUrl=queueUrl, MessageBody=message)
#    end

    @showprogress 1 "producing tasks..." for i in 1:10:length(messageList)
        messages = messageList[i:min(i+9, length(messageList))]
        messageBatch = map((x,y)->["Id"=>string(i+x-1), "MessageBody"=>y],
                                                1:length(messages), messages)
        SQS.send_message_batch(AWS_CREDENTIAL; QueueUrl=queueUrl, 
                               SendMessageBatchRequestEntry=messageBatch)
    end 

    nothing
end

main()
