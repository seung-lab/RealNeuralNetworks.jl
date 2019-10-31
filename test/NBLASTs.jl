using Test 

using RealNeuralNetworks.Neurons 
using RealNeuralNetworks.NBLASTs
using RealNeuralNetworks.Utils.VectorClouds 
using RealNeuralNetworks.Utils.RangeIndexingArrays 
using RealNeuralNetworks.Neurons.Segments 

using CSV

ASSET_DIR = joinpath(@__DIR__, "../asset/")

NEURON_ID1 = 77625
NEURON_ID2 = 77641

const k = 20

@testset "test NBLAST module..." begin 
    neuron1 = Neurons.load(joinpath(ASSET_DIR, "$(NEURON_ID1).swc.bin"));
    neuron2 = Neurons.load(joinpath(ASSET_DIR, "$(NEURON_ID2).swc.bin"));


    vectorCloud1 = NBLASTs.VectorCloud(neuron1)
    vectorCloud2 = NBLASTs.VectorCloud(neuron2)

    # transform to micron
    vectorCloud1[1:3, :] ./= 1000
    vectorCloud2[1:3, :] ./= 1000

    # read the precomputed score matrix, which is the joint distribution of the scores
    println("\nread the precomputed score table...")
    SMAT_FCWB_PATH = joinpath(ASSET_DIR, "smat_fcwb.csv")
    df = CSV.read(SMAT_FCWB_PATH)    
    ria = RangeIndexingArray{Float32}(df)

    println("\ncompute nblast score...")
    @time score = NBLASTs.nblast(vectorCloud2, vectorCloud1; ria=ria) 
    @time score = NBLASTs.nblast(vectorCloud2, vectorCloud1; ria=ria) 
    println("query $(NEURON_ID1) against target $(NEURON_ID2): ", score)
    @test isapprox(score, Float32(50887.82f0))
    
    println("\ncompute nblast score...")
    @time score = NBLASTs.nblast(vectorCloud1, vectorCloud2; ria=ria) 
    println("query $(NEURON_ID1) against target $(NEURON_ID2): ", score)
    @test isapprox(score, Float32(53697.95f0))

    #@time score = NBLASTs.nblast(vectorCloud1, vectorCloud1; ria=ria)
    #println("query $(NEURON_ID1) against itself: ", score)
    #@test isapprox(score, Float32(86507.24f0)) 


    neuronList = [neuron1, neuron2]
    # the result from R NBLAST is :
    # ID    77625	    77641
    # 77625	86501.20 	53696.72
    # 77641	50891.03 	101011.08
    RNBLAST_RESULT = Float32[86501.20 	53696.72; 50891.03 	101011.08]
    println("\ncompute similarity matrix...")
    @time rawSimilarityMatrix, normalizedSimilarityMatrix, meanSimilarityMatrix = NBLASTs.nblast_allbyall(neuronList; ria=ria, k=k,
                                                        downscaleFactor=1000)
    @show rawSimilarityMatrix
    @test isapprox.(rawSimilarityMatrix,  RNBLAST_RESULT) |> all
    

    #println("\ncompute normalised similarity matrix...")
    #RNBLAST_RESULT = Float32[1.0 0.5315924; 0.5883275 1.0]
    #@show normalizedSimilarityMatrix
    #@test isapprox.(normalizedSimilarityMatrix,  RNBLAST_RESULT) |> all
    #
    #RNBLAST_RESULT = Float32[1.0 0.5599599; 0.5599599 1.0]
    #println("\ncompute mean similarity matrix...")
    #@show meanSimilarityMatrix
    #@test isapprox.(meanSimilarityMatrix,  RNBLAST_RESULT) |> all
    #
    #vectorCloudList = [vectorCloud1, vectorCloud2]
    #RNBLAST_RESULT = Float32[86501.20 	53696.72; 50891.03 	101011.08]
    #@time rawSimilarityMatrix = NBLASTs.nblast_allbyall(vectorCloudList; ria=ria)
    #@show rawSimilarityMatrix
    #@test isapprox.(rawSimilarityMatrix,  RNBLAST_RESULT) |> all
 
    #    
    #println("\ncompute nblast score using default zfish score table ...")
    ## this test will not work since neuron reading from swc/swc.bin 
    ## do not support node class yet
    #vectorCloud1 = NBLASTs.VectorCloud(neuron1; class=Segments.DENDRITE_CLASS)  
    #vectorCloud2 = NBLASTs.VectorCloud(neuron2; class=Segments.DENDRITE_CLASS)  
    #@time score = NBLASTs.nblast(vectorCloud2, vectorCloud1;) 
    #println("query $(NEURON_ID1) against target $(NEURON_ID2): ", score)
    ## @test isapprox(score, -63722.48f0)
    #
    #println("\ntest similarity measurement for dendrites...")
    #@time rawSimilarityMatrix, normalizedSimilarityMatrix, meanSimilarityMatrix = 
    #            NBLASTs.nblast_allbyall(neuronList; ria=ria, k=k, downscaleFactor=1000)
    #@show rawSimilarityMatrix
    ##println("similarity matrix of dendrites: ", similarityMatrix)

    #println("compute similarity score with semantic small-to-big...")
    #semanticSmall2bigSimilarityMatrix = 
    #            NBLASTs.nblast_allbyall_small2big(neuronList; k=k, ria=ria, semantic=true);

end # end of test module
