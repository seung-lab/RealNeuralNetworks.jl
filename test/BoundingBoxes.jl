using Base.Test
using RealNeuralNetworks.BranchNets.Branches.BoundingBoxes 

function create_fake_bounding_box()
    bbox = BoundingBox((1,1,1), (204,204,300))
end 

@testset "test bounding box" begin
    bbox = BoundingBox((1,2,3), (2,3,4))
    d = BoundingBoxes.distance_from(bbox, (3,4,5))
    @test_approx_eq(d, sqrt(3))

    bbox2 = BoundingBox((2,3,4), (3,4,5))
    bbox3 = union(bbox, bbox2)
    @show BoundingBox((1,2,3), (3,4,5))
    @show bbox3
    @test isequal(bbox3, BoundingBox((1,2,3), (3,4,5)))
    @test bbox3 == BoundingBox((1,2,3), (3,4,5))
end 
