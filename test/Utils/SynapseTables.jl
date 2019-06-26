using Test
using CSV
using DataFrames
using RealNeuralNetworks.Utils.SynapseTables 

@testset "synapse tables" begin
    syn = CSV.read(joinpath(@__DIR__, "../../asset/syn100.csv"); copycols=true)
    @test 100 == DataFrames.nrow(syn)
    
    println("\npreprocessing...")
    @time syn = SynapseTables.preprocess(syn, (5,5,45))
    @test 100 >= DataFrames.nrow(syn)

    syn1 = SynapseTables.get_synapses_of_a_neuron( syn, 76263)
    syn2 = SynapseTables.get_pre_synapses_of_a_neuron(syn, 76782)
    syn3 = SynapseTables.get_post_synapses_of_a_neuron(syn, 76263)
    @test DataFrames.nrow(syn1) > 0
    @test DataFrames.nrow(syn2) > 0
    @test DataFrames.nrow(syn3) > 0

    bbox = SynapseTables.BoundingBox( syn, "psd_" )

    presynCoordinates = SynapseTables.get_coordinate_array(syn, "presyn_")
    @test !isempty( presynCoordinates )

    mask = SynapseTables.get_mask( syn )
    @test any(mask .> zero(eltype(mask)) )

    println("\npostprocessing...")
    @time SynapseTables.postprocess(syn, (5,5,45))
end 
