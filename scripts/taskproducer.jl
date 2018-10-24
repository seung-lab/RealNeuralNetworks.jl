#!/usr/bin/env julia
using JSON
using AWSSDK.SQS 
using AWSCore
using ProgressMeter
using ArgParse 

const AWS_CREDENTIAL = AWSCore.aws_config()


function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--jsonfile", "-j"
            help = "the id list was contained in the consensus json file"
            arg_type = String
        "--sqsqueue", "-q"
            help = "AWS SQS queue name"
            arg_type = String 
            default = "zfish"
    end 
    return parse_args(s)
end 


function main()
    
    args = parse_commandline()
    @show args
     
    d = JSON.parsefile(args["jsonfile"] |> expanduser) 
    neuronIdSet = Set{Int}()
    for v in d |> values  
        for neuronIdStr in keys(v) 
            neuronId = Meta.parse(neuronIdStr)
            push!(neuronIdSet, neuronId)
        end
    end 

    messageList = Vector{String}()
    for id in neuronIdSet 
        push!(messageList, string(id))
    end

    # put to AWS SQS queue
    queueUrl = SQS.get_queue_url( QueueName=args["sqsqueue"] )["QueueUrl"]

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
