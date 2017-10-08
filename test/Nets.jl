using Base.Test
using TEASAR.Nets
using TEASAR.Skeletons

const ONE_UINT32 = UInt32(1)

function create_cylinder_segmentation(radius::Float32=Float32(20), height::Integer=50)
    width = Int( radius*2 + 1 )
    center = Int(radius + 1)
    seg = zeros(UInt32, (width,width,height))
    for x in Int(center-radius):Int(center+radius)
        for y in Int(center-radius):Int(center+radius)
            if norm([x,y].-[center, center]) < radius 
                seg[x,y,:] = ONE_UINT32 
            end 
        end 
    end 
    seg
end

@testset "test Nets" begin 
    seg = create_cylinder_segmentation()
    println("skeletonization ...")
    @time skeleton = Skeleton(seg)
    println("transform to Net structure")
    @time net = Net(skeleton)
end 
