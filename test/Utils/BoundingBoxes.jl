using Test
using RealNeuralNetworks.Utils.BoundingBoxes 

function create_fake_bounding_box()
    bbox = BoundingBox((1,1,1), (204,204,300))
end 

@testset "test bounding box" begin
    bbox = BoundingBox((1,2,3), (2,3,4))
    d = BoundingBoxes.distance_from(bbox, (3,4,5))
    @test d ≈ sqrt(3)
    range = BoundingBoxes.get_unit_range(bbox)
    @test range == (1:2, 2:3, 3:4)
    @test size(bbox) == (2,2,2)


    bbox2 = BoundingBox((2,3,4), (3,4,5))
    bbox3 = union(bbox, bbox2)
    @show BoundingBox((1,2,3), (3,4,5))
    @show bbox3
    @test isequal(bbox3, BoundingBox((1,2,3), (3,4,5)))
    @test bbox3 == BoundingBox((1,2,3), (3,4,5))
end 
