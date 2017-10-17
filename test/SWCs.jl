using RealNeuralNetworks.SWCs
using Base.Test

@testset "test SWC" begin 
    # read swc
    exampleFile = joinpath(dirname(@__FILE__), "../assert/example.swc")
    swc = SWCs.load( exampleFile )
    str = String(swc)
    tempFile = tempname() * ".swc"
    SWCs.save(swc, tempFile)
    @test readstring(exampleFile) == readstring( tempFile )
    rm(tempFile)

    # test stretch
    # println("test stretch according to voxel size...")
    SWCs.stretch_coordinates!(swc, (2,3,4))
end 

