using RealNeuralNetworks
using RealNeuralNetworks.NodeNets.PointArrays

using Test

@testset "test PointArray add_offset! operation" begin
    # build fake data
    a = [   0xc95f39a7  0x21a2a09c  0xc9e447a2
            0xa8e1555e  0x7fbeb45b  0xc786dc67
            0x940649cd  0x20b148bc  0xc15152cb
            0xc91b9db1  0x39fee982  0x7019e841]
    offset = (UInt32(1), UInt32(2), UInt32(3))
    PointArrays.add_offset!(a, offset)
    b = [0xc95f39a8  0x21a2a09e  0xc9e447a5
         0xa8e1555f  0x7fbeb45d  0xc786dc6a
         0x940649ce  0x20b148be  0xc15152ce
         0xc91b9db2  0x39fee984  0x7019e844]
    @test all(a.==b)
end 
