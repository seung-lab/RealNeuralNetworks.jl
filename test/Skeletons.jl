using TEASAR  
using Base.Test
using HDF5
using BigArrays
using GSDicts
using TEASAR.Skeletons
using TEASAR.SWCs

const DEFAULT_VOXEL_SIZE = (UInt32(80), UInt32(80), UInt32(40))

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

function get_seg_from_gs()
    ba = BigArray(GSDict("gs://neuroglancer/zfish_v1/consensus-20170829/80_80_45"))
    seg = ba[2457:2968, 1777:2288, 16401:16912]
end 

function create_fake_seg()
    seg = zeros(UInt32,(100,100,100))
    seg[50,50,:] = UInt32(1)
    seg[49:52, 49:52, 48:52] = 1
    seg[47:54, 47:54, 71:78] = 1
    return seg 
end 

@testset "test skeletonization" begin
    seg = create_fake_seg()
    # @time seg = get_seg_from_h5()
    # @time seg = get_seg_from_gs()
    # @time swc = TEASAR.skeletonize(seg; voxel_size=DEFAULT_VOXEL_SIZE)
    println("building skeleton ...")
    # @time skeleton = Skeleton( seg; obj_id = UInt32(76880) )
    @time skeleton = Skeleton( seg; obj_id = UInt32(1) )
    bin = Skeletons.get_neuroglancer_precomputed(skeleton)
    open("/tmp/fake.bin", "w") do f
        write(f, bin)
    end 
    @time swc = SWC( skeleton )
    @test length(swc) > 1
    SWCs.save(swc, tempname() * ".swc")
end 


