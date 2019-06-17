#!/usr/bin/env julia
using ProgressMeter
using ArgParse

using RealNeuralNetworks
using RealNeuralNetworks.SWCs
using RealNeuralNetworks.Neurons 

function parse_commandline()
    s = ArgParseSettings()
    @add_arg_table s begin 
        "--inputdir", "-i"
            help = "the input directory to store swc files"
            arg_type = String
            default = "/tmp"
        "--outputpath", "-o"
            help = "the output directory to store binary swc files"
            arg_type = String
            default = "/tmp"
        "--resample-node-distance", "-d"
            help = "mip level of the dataset" 
            arg_type = Float32
            default = Float32(100)
    end
    return parse_args(s)
end 
 
function main()
    args = parse_commandline()

    @showprogress 1 "postprocessing..." for fileName in readdir(args["inputdir"])
        fullName = joinpath(expanduser(args["inputdir"]), fileName)
        if endswith(fileName, "swc.bin")
            neuron = Neurons.load_swc_bin(fullName)
        elseif endswith(fileName, ".swc")
            neuron = Neurons.load_swc(fullName)
        elseif isdir(fullName)
            warn("this is a directory: $(fileName)")
        else 
            warn("unsupported file format: $(fileName)")
            continue 
        end 
        neuron = Neurons.remove_hair(neuron)
        neuron = Neurons.remove_redundent_nodes(neuron)
        neuron = Neurons.remove_subtree_in_soma(neuron)
        neuron = Neurons.remove_terminal_blobs(neuron)
        neuron = Neurons.resample(neuron, args["resample-node-distance"])

        # output file 
        outputName = joinpath(expanduser(args["outputpath"]), fileName)
        if !isdir(dirname(outputName))
            mkdir(dirname(outputName))
        end 
        Neurons.save_swc_bin(neuron, outputName)
    end
end

main()
