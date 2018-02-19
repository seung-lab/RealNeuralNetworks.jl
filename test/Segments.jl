using RealNeuralNetworks
using RealNeuralNetworks.Neurons
using RealNeuralNetworks.Neurons.Segments
using Base.Test

const SWC_BIN_PATH = joinpath(@__DIR__, "../assert/77625.swc.bin")

@testset "test Segments" begin
    # construct a segment
    neuron = Neurons.load_swc_bin( SWC_BIN_PATH )
    # get a random segment
    segment = neuron[5]
    
    println("get tortuosity...")
    @show Segments.get_tortuosity( segment )
    println("get tail head radius ratio ...")
    @show Segments.get_tail_head_radius_ratio( segment )
end 
