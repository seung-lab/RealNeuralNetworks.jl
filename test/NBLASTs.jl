using Test 

using RealNeuralNetworks.SWCs
using RealNeuralNetworks.NodeNets 
using RealNeuralNetworks.Neurons 
using RealNeuralNetworks.NBLASTs
using RealNeuralNetworks.Utils.VectorClouds 
using RealNeuralNetworks.Utils.RangeIndexingArrays 
using RealNeuralNetworks.Neurons.Segments 

using CSV

ASSET_DIR = joinpath(@__DIR__, "../asset/")

NEURON_ID1 = 77625
NEURON_ID2 = 77641


function prepare_neurons()
    neuron1 = Neurons.load("../01_data/0920/postprocessed/atlas_space/bin/$(NEURON_ID1).bin")
    neuron2 = Neurons.load("../01_data/0920/postprocessed/atlas_space/bin/$(NEURON_ID2).bin")

    Neurons.save(neuron1, "/tmp/$(NEURON_ID1).swc.bin")
    Neurons.save(neuron2, "/tmp/$(NEURON_ID2).swc.bin")
end 

@testset "test NBLAST module..." begin 
    neuron1 = SWCs.load(joinpath(ASSET_DIR, "$(NEURON_ID1).swc.bin")) |> NodeNet;
    neuron2 = SWCs.load(joinpath(ASSET_DIR, "$(NEURON_ID2).swc.bin")) |> NodeNet;

    k = 20

    vectorCloud1 = NBLASTs.VectorCloud(neuron1)
    vectorCloud2 = NBLASTs.VectorCloud(neuron2)

    println("\ntest data structure transformation...")
    NBLASTs.VectorCloud(Neuron(neuron1))
    NBLASTs.VectorCloud(neuron1)

    # transform to micron
    vectorCloud1[1:3, :] ./= 1000
    vectorCloud2[1:3, :] ./= 1000

    # read the precomputed score matrix, which is the joint distribution of the scores
    println("\nread the precomputed score table...")
    SMAT_FCWB_PATH = joinpath(@__DIR__, "../asset/smat_fcwb.csv")
    df = CSV.read(SMAT_FCWB_PATH)    
    ria = RangeIndexingArray{Float32}(df)

    println("\ncompute nblast score...")
    @time score = NBLASTs.nblast(vectorCloud2, vectorCloud1; ria=ria) 
    println("query $(NEURON_ID1) against target $(NEURON_ID2): ", score)
    @test isapprox(score, Float32(43072.21))
    @time score = NBLASTs.nblast(vectorCloud1, vectorCloud1; ria=ria)
    println("query $(NEURON_ID1) against itself: ", score)
    @test isapprox(score, Float32(69314.85)) 


    vectorCloudList = [vectorCloud1, vectorCloud2]
    # the result from R NBLAST is :
    # ID    77625	    77641
    # 77625	69314.85	45313.29
    # 77641	43072.21	80601.58
    RNBLAST_RESULT = Float32[69314.85 45313.29; 43072.21 80601.58]
    println("\ncompute similarity matrix...")
    @time similarityMatrix = NBLASTs.nblast_allbyall(vectorCloudList; ria=ria, normalisation=:raw)
    @show similarityMatrix
    @test isapprox.(similarityMatrix,  RNBLAST_RESULT) |> all
 
    println("\ncompute normalised similarity matrix...")
    RNBLAST_RESULT = Float32[1.0 0.5621886; 0.6213994 1.0]
    @time similarityMatrix = NBLASTs.nblast_allbyall(vectorCloudList; ria=ria, normalisation=:normalised)
    @show similarityMatrix
    @test isapprox.(similarityMatrix,  RNBLAST_RESULT) |> all

    
    RNBLAST_RESULT = Float32[1.0 0.591794; 0.591794 1.0]
    println("\ncompute mean similarity matrix...")
    @time similarityMatrix = NBLASTs.nblast_allbyall(vectorCloudList; ria=ria, normalisation=:mean)
    @show similarityMatrix
    @test isapprox.(similarityMatrix,  RNBLAST_RESULT) |> all

    println("\ntest similarity measurement for dendrites...")
    vectorCloud1 = NBLASTs.VectorCloud(neuron1; class=Segments.DENDRITE_CLASS)  
    vectorCloud2 = NBLASTs.VectorCloud(neuron2; class=Segments.DENDRITE_CLASS)  
    vectorCloudList = [vectorCloud1, vectorCloud2] 
    @time similarityMatrix = NBLASTs.nblast_allbyall(vectorCloudList; ria=ria, normalisation=:raw)
    @show similarityMatrix
    #println("similarity matrix of dendrites: ", similarityMatrix)
    
    println("\ncompute nblast score using default zfish score table ...")
    @time score = NBLASTs.nblast(vectorCloud2, vectorCloud1;) 
    println("query $(NEURON_ID1) against target $(NEURON_ID2): ", score)
    @test isapprox(score, -51068.773) 
end 
