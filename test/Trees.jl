using Base.Test
using TEASAR.Skeletons 
using TEASAR.Trees

function create_cylinder_seg(sz::Tuple=(51,51,101), radius::Float32=Float32(30))
    center = map(x->div(x,2), sz)
    seg = zeros(UInt32, sz)
    for y in 1:sz[2]
        for x in 1:sz[1]
            if norm( [x,y] .- [center[1:2]...] ) < radius 
                seg[x,y,:] = UInt32(1)
            end 
        end 
    end 
    seg
end 

@testset "test Trees" begin
    seg = create_cylinder_seg()
    skeleton = Skeleton(seg)
    tree = Tree(skeleton)
end # testset 
