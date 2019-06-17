using RealNeuralNetworks
using RealNeuralNetworks.Neurons
using RealNeuralNetworks.Neurons.Segments
using Test

const SWC_BIN_PATH = joinpath(@__DIR__, "../asset/77625.swc.bin")

@testset "test Segments" begin
    # construct a segment
    neuron = Neurons.load_swc_bin( SWC_BIN_PATH )
    # get a random segment
    println("indexing from neuron...")
    segment = neuron[5]
    
    println("get tortuosity...")
    @show Segments.get_tortuosity( segment )
    println("get tail head radius ratio ...")
    @show Segments.get_tail_head_radius_ratio( segment )

    println("merge segments...")
    segment2 = neuron[6]
    merged_segment = merge(segment, segment2)
    @test length(segment) + length(segment2) == length(merged_segment)
    
    println("remove some nodes...")
    newSegment = Segments.remove_nodes(segment, 2:4)
    @test length(newSegment) == length(segment) - 3

    println("remove redundent nodes...")
    Segments.remove_redundent_nodes!(segment) 
end 
