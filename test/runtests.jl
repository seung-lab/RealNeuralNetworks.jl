using TEASAR.Skeleton  
using Base.Test
using HDF5

function get_seg_from_h5()
    # read seg data
    f = h5open("/usr/people/jingpeng/seungmount/research/Ashwin/Scripts/NG_scripts/77605.h5")
    seg = f["main"][1,:,:,101:500]
    close(f)
    seg = reshape(seg, size(seg)[2:4])
    @assert ndims(seg) == 3
    for i in eachindex(seg)
        if seg[i] == 77605
            seg[i] = convert(UInt32,1)
        else
            seg[i] = convert(UInt32, 0)
        end 
    end 
    return seg 
end

seg = zeros(UInt32,(100,100,100))
seg[50,50,:] = UInt32(1)
seg[49:52, 49:52, 48:52] = 1
seg[47:54, 47:54, 71:78] = 1

#@time seg = get_seg_from_h5()

@testset "test teasar" begin 
    @time points, edges, nodes, roots, radii, dests = Skeleton.skeletonize(seg)
    @test !isempty(edges)
    @show nodes
    @show edges 
    @show roots 
    @show radii 
    @show dests
end 
