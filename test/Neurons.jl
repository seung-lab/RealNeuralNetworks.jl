using RealNeuralNetworks  
using Test


@testset "test Neuron type methods..." begin 
    # read swc
    exampleFile = joinpath(@__DIR__, "../asset/Nov10IR3e.CNG.swc")
    println("load plain text swc...")
    neuron = load_swc( exampleFile )
    @time neuron = load_swc( exampleFile )
    println("loaded swc with nodes: ", length(neuron))

    println("compute total path length")
    pathLength = get_path_length(neuron)
    println("path length of this neuron: ", pathLength)
    # we need to verify the L-measure number using other tools
    # matlab treestoolbox gives 3516.6
    # L-measure result: 3505.76
    @test  isapprox(pathLength, 3516.6; atol = 0.1)

    tempFile = tempname() * ".swc"
    println("save plain text swc ...")
    @time save_swc(neuron, tempFile)
    rm(tempFile)
    println("save binary neuron...")
    @time save_neuron_bin( neuron, "$(tempFile).neuron.bin" )
    println("load binary neuron...")
    @time neuron2 = load_neuron_bin( "$(tempFile).neuron.bin" )
    @test neuron == neuron2
    rm("$(tempFile).neuron.bin")
 
    println("add offset...")
    add_offset!(neuron, (-one(Float32),-one(Float32),-one(Float32)))
    println("get neuroglancer precomputed...")
    bin = get_neuroglancer_precomputed(neuron)
    # open("/tmp/fake.bin", "w") do f write(f, bin)  end 
    println("stretch coordinates...")
    stretch_coordinates!(neuron, 4)
    @test length(neuron) > 1
    
    # test stretch
    stretch_coordinates!(neuron, (2,3,4))
end 
