using RealNeuralNetworks  
using Test


@testset "test NodeNet type methods..." begin 
    # read swc
    exampleFile = joinpath(@__DIR__, "../asset/Nov10IR3e.CNG.swc")
    println("load plain text swc...")
    @time nodeNet = load_swc( exampleFile )

    tempFile = tempname() * ".swc"

    println("save plain text swc ...")
    @time save_swc(nodeNet, tempFile)
    #@test read(exampleFile, String) == read( tempFile , String)
    rm(tempFile)
    
    println("save binary nodenet ...")
    @time save_nodenet_bin( nodeNet, "$(tempFile).nodenet.bin" )
    println("load binary nodenet ...")
    @time nodeNet2 = load_nodenet_bin( "$(tempFile).nodenet.bin" )
    @test nodeNet == nodeNet2
    rm("$(tempFile).nodenet.bin")
 
    println("add offset...")
    add_offset!(nodeNet, (-one(Float32),-one(Float32),-one(Float32)))
    println("get neuroglancer precomputed...")
    bin = get_neuroglancer_precomputed(nodeNet)
    # open("/tmp/fake.bin", "w") do f write(f, bin)  end 
    println("stretch coordinates...")
    stretch_coordinates!(nodeNet, 4)
    @test length(nodeNet) > 1
    
    println("compute total path length")
    pathLength = get_path_length(nodeNet)
    println("path length of this neuron: ", pathLength)
    # we need to verify the L-measure number using other tools
    # @test  isapprox(pathLength, 3505.76; atol = 0.1)

    # test stretch
    stretch_coordinates!(nodeNet, (2,3,4))
end 