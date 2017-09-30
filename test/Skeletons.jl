using TEASAR  
using Base.Test
using HDF5
using BigArrays
using GSDicts
using TEASAR.Skeletons
using TEASAR.SWCs
using JLD

const CELL_ID = UInt32(76880)
const OFFSET = (UInt32(2456), UInt32(1776), UInt32(16400))
const EXPANSION= (UInt32(80), UInt32(80), UInt32(40))
const GS_SEG_PATH = "gs://neuroglancer/zfish_v1/consensus-20170829/80_80_45"
const GS_SKELETON_PATH = "gs://neuroglancer/zfish_v1/consensus-20170829/skeleton_mip_4"

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
    ba = BigArray(GSDict( GS_SEG_PATH ))
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
    println("building skeleton ...")
    #@time skeleton = Skeleton( seg; obj_id = CELL_ID )
    # Skeletons.add_offset!(skeleton, OFFSET)
    @time skeleton = Skeleton( seg; obj_id = UInt32(1) )
    bin = Skeletons.get_neuroglancer_precomputed(skeleton)
    # open("/tmp/fake.bin", "w") do f write(f, bin)  end 
    @time swc = SWC( skeleton )
    @test length(swc) > 1
    save("/tmp/$(CELL_ID).jld", "skeleton", skeleton, "swc", swc)
    SWCs.save(swc, "/tmp/$(CELL_ID).swc")

    # test saving to google cloud for neuroglancer visualization
    d_json = GSDict( GS_SKELETON_PATH; valueType = String )
    d_bin  = GSDict( GS_SKELETON_PATH )
    Skeletons.save(skeleton, CELL_ID, d_json, d_bin)
end 


