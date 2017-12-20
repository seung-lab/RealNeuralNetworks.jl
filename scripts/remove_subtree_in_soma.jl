#!/usr/bin/env julia

using RealNeuralNetworks
using RealNeuralNetworks.SWCs
using RealNeuralNetworks.Neurons 

include("Common.jl"); using Common

const CELL_ID = 77390

function main()
    args = parse_commandline()
    swc = SWCs.load_gzip_swc("swcgz/$(CELL_ID).swc.gz")
    neuron = Neuron( swc )
    neuron = Neurons.remove_subtree_in_soma(neuron)
    neuron = Neurons.remove_hair( neuron )
    Neurons.save(neuron, "$(CELL_ID).swc")
end

main()
