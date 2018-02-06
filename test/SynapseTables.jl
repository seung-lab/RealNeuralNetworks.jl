using Base.Test
using CSV
using DataFrames
using RealNeuralNetworks.SynapseTables 

@testset "synapse tables" begin
    syn = CSV.read(joinpath(@__DIR__, "../assert/syn100.csv"))
    @test 100 == DataFrames.nrow(syn)
    SynapseTables.preprocessing!(syn)
    @test 100 >= DataFrames.nrow(syn)

    syn1 = SynapseTables.get_synapses_of_a_neuron( syn, 76263)
    syn2 = SynapseTables.get_synaptic_boutons_of_a_neuron(syn, 76782)
    syn3 = SynapseTables.get_postsynaptic_density_of_a_neuron(syn, 76263)
    @test DataFrames.nrow(syn1) > 0
    @test DataFrames.nrow(syn2) > 0
    @test DataFrames.nrow(syn3) > 0

    bbox = SynapseTables.BoundingBox( syn )
end 
