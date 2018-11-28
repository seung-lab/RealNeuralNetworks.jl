using Test 

using RealNeuralNetworks.NeuralNets 
using RealNeuralNetworks.Neurons

using CSV

@testset "test NeuralNets module..." begin 
    neuronId2neuron = Dict{Int, Neuron}()
    neuronId2neuron[77625] = Neurons.load(joinpath(__DIR__, "../asset/77625.swc.bin"))
    neuronId2neuron[77641] = Neurons.load(joinpath(__DIR__, "../asset/77641.swc.bin"))

    df = CSV.load(joinpath(__DIR__, "../asset/77265.pre.synapses.csv")
    #Neurons.attach_pre_synapses()
    net = NeuralNet(neuronId2neuron)
end 
