using Test 
using CSV 
using DataFrames 
using RealNeuralNetworks.Neurons.Segments.Synapses 

@testset "test Synapses module" begin 
    df = CSV.read(joinpath(Base.@__DIR__, 
                           "../asset/77625.pre.synapses.csv"))
    synapseList = Vector{Synapse{Float32}}()
    for row in DataFrames.eachrow(df)
        push!(synapseList, Synapse(row))
    end 
    @test !isempty( synapseList ) 

    synapse = synapseList[1]
    @test Synapses.get_psd_size(synapse) > 0 
    @test Synapses.get_psd_segmentation_id(synapse) > 0
    @test Synapses.get_psd_coordinate( synapse ) |> length > 0
    Synapses.get_psd_bounding_box( synapse ) 

    @test Synapses.get_pre_synaptic_segmentation_id( synapse ) > 0
    @test Synapses.get_pre_synaptic_coordinate( synapse ) |> length > 0
    @test Synapses.get_pre_synaptic_weight( synapse ) > 0.0 

    @test Synapses.get_post_synaptic_segmentation_id( synapse ) > 0
    @test Synapses.get_post_synaptic_coordinate( synapse ) |> length > 0
    @test Synapses.get_post_synaptic_weight( synapse )  > 0.0
end 
