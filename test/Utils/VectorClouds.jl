using Test 
using RealNeuralNetworks.SWCs
using RealNeuralNetworks.NodeNets 
using RealNeuralNetworks.NBLASTs 
using RealNeuralNetworks.Utils.VectorClouds 

ASSET_DIR = joinpath(@__DIR__, "../../asset/")

NEURON_ID1 = 77641
NEURON_ID2 = 77625

@testset "test NBLAST module..." begin 
    neuron1 = SWCs.load(joinpath(ASSET_DIR, "$(NEURON_ID1).swc.bin")) |> NodeNet;
    neuron2 = SWCs.load(joinpath(ASSET_DIR, "$(NEURON_ID2).swc.bin")) |> NodeNet; 
    
    vectorCloud1 = NBLASTs.VectorCloud(neuron1)
    vectorCloud2 = NBLASTs.VectorCloud(neuron2)
    
    distMatrix = VectorClouds.query(vectorCloud1, vectorCloud2)
    @test any(distMatrix .!== zero(Float32))
end 
