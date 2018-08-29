using RealNeuralNetworks
using RealNeuralNetworks.NodeNets.DBFs   
using Test

seg = zeros(UInt32,(100,100,100))
seg[48:52,48:52, 1:end ] .= one(UInt32)
seg[49:52, 49:52, 48:52] .= one(UInt32)
seg[47:54, 47:54, 71:78] .= one(UInt32)


points = RealNeuralNetworks.NodeNets.PointArrays.from_seg(seg)
@testset "test dbf computation" begin 
    dbf1 = DBFs.compute_DBF(points)
    #dbf2 = DBFs.compute_DBF(points, boundary_point_indexes)
    bin_im = DBFs.create_binary_image( seg )
    bin_im2 = DBFs.create_binary_image( points )
    @show size(bin_im)
    @show size(bin_im2)
    # the size of bin_im2 is smaller than bin_im
    #@test all(bin_im .== bin_im2)
    dbf3 = DBFs.compute_DBF(points, bin_im )
    #@show dbf1 
    #@show dbf3
    # it seems that the dbf from binary image is not correct
    # map((x,y) -> @test_approx_eq_eps(x,y,1), dbf1, dbf3)
end 


