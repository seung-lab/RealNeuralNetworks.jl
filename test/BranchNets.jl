using Base.Test
using TEASAR.BranchNets
using TEASAR.NodeNets
using TEASAR.SWCs

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

@testset "test BranchNets" begin 
    seg = create_cylinder_segmentation()
    println("skeletonization to build a NodeNet ...")
    @time nodeNet = NodeNet(seg)
    @show nodeNet
    println("transform to BranchNet structure ...")
    @time branchNet = BranchNet(nodeNet)
    @show branchNet 
    println("transform to SWC structure ...")
    @time swc = SWCs.SWC( branchNet )
    @show swc
    SWCs.save(swc, "/tmp/cylinder.swc")
end 
