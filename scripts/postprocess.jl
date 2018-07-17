#!/usr/bin/env julia

using RealNeuralNetworks
using RealNeuralNetworks.SWCs
using RealNeuralNetworks.Neurons 

include("Common.jl"); using Common

function main()
    args = parse_commandline()

    for filename in dirname(args["swcdir"])
        @show filename
        neuron = Neurons.load_swc_bin(joinpath(args["swcdir"], filename))
        neuron = Neurons.remove_hair(neuron)
        neuron = Neurons.remove_redundent_nodes(neuron)
        neuron = Neurons.remove_subtree_in_soma(neuron)
        neuron = Neurons.remove_terminal_blobs(neuron)
        neuron = Neurons.resample(neuron, Float32(100))
        Neurons.save_swc_bin(neuron, joinpath(args["outputpath"], filename))
        Neurons.save_swc(neuron, joinpath(args["outputpath"]*"_swc", filename))
    end
end

main()
