using TEASAR
using Base.Test

seg = zeros(UInt32,(100,100,100))
seg[48:52,48:52,:] = UInt32(1)
seg[49:52, 49:52, 48:52] = 1
seg[47:54, 47:54, 71:78] = 1

points, boundary_point_indexes = TEASAR.PointArrays.from_seg(seg)

@testset "test dbf computation" begin 
    dbf1 = TEASAR.DBFs.compute_DBF(points)
    dbf2 = TEASAR.PointArrays.compute_DBF(points, boundary_point_indexes)
    @show dbf1 
    @show dbf2
    map((x,y) -> @test_approx_eq_eps(x,y,1), dbf1, dbf2)
end 

