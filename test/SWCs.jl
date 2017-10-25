using RealNeuralNetworks.SWCs
using Base.Test
using Libz

@testset "test SWC" begin 
    # read swc
    exampleFile = joinpath(dirname(@__FILE__), "../assert/example.swc")
    println("load plain text swc ...")
    @time swc = SWCs.load( exampleFile )
    str = String(swc)
    tempFile = tempname() * ".swc"

    println("save plain text swc ...")
    @time SWCs.save(swc, tempFile)
    @test readstring(exampleFile) == readstring( tempFile )
    rm(tempFile)

    println("save gzip swc file: $(tempFile).gz ...")
    @time SWCs.save_gzip_swc(swc, "$(tempFile).gz")
    stream = open( "$(tempFile).gz" ) |> ZlibInflateInputStream
    str2 = readstring( stream )
    @test str == str2

    println("load gzip swc ...")
    @time swc2 = SWCs.load_gzip_swc( "$(tempFile).gz" )
    @test length(swc) == length(swc2)
    @test String(swc) == String(swc2)
    # @test swc == swc2 
    rm("$(tempFile).gz")
    
    # test stretch
    SWCs.stretch_coordinates!(swc, (2,3,4))
end 

