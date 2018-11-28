using CSV 
using RealNeuralNetworks 
using RealNeuralNetworks.Utils.RangeIndexingArrays 
using Test 


@testset "test RangeIndexingArrays module..." begin 
    SMAT_FCWB_PATH = joinpath(@__DIR__, "../../asset/smat_fcwb.csv")
    df = CSV.read(SMAT_FCWB_PATH)    
    
    ria = RangeIndexingArray{Float32}(df)
    @test ria[Float32(0.2), Float32(0.03)] == 9.500097f0
    @test ria[Float32(9.2), Float32(0.25)] == 1.1505605f0
end 
