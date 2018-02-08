using Base.Test
using CSV
using DataFrames
using RealNeuralNetworks.SynapseTables 

@testset "synapse tables" begin
    syn = CSV.read(joinpath(@__DIR__, "../assert/syn100.csv"))
    @test 100 == DataFrames.nrow(syn)
    SynapseTables.preprocessing!(syn, (5,5,45))
    @test 100 >= DataFrames.nrow(syn)

    syn1 = SynapseTables.get_synapses_of_a_neuron( syn, 76263)
    syn2 = SynapseTables.get_synaptic_boutons_of_a_neuron(syn, 76782)
    syn3 = SynapseTables.get_postsynaptic_density_of_a_neuron(syn, 76263)
    @test DataFrames.nrow(syn1) > 0
    @test DataFrames.nrow(syn2) > 0
    @test DataFrames.nrow(syn3) > 0

    bbox = SynapseTables.BoundingBox( syn, "COM" )

    presynCoordinates = SynapseTables.get_coordinate_array(syn, "presyn")
    @test !isempty( presynCoordinates )

    mask = SynapseTables.get_mask( syn )
    @test any(mask .> zero(eltype(mask)) )
end 
