using Test

using RealNeuralNetworks  
using BigArrays
using BigArrays.GSDicts 

using LinearAlgebra
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

function get_seg_from_gs()
    ba = BigArray(GSDict( SEG_PATH ))
    seg = ba[2457:2968, 1777:2288, 16401:16912]
end 


function create_fake_seg()
    seg = zeros(UInt32,(100,100,100))
    seg[50,50,:] .= one(UInt32)
    seg[49:52, 49:52, 48:52] .= one(UInt32)
    seg[47:54, 47:54, 71:78] .= one(UInt32)
    return seg 
end 

@testset "test TEASAR" begin
    seg = create_fake_seg()
    # @time seg = get_seg_from_h5()
    # @time seg = get_seg_from_gs()
    println("\nbuilding nodeNet ...")
    #@time nodeNet = NodeNet( seg; obj_id = CELL_ID )
    @time nodeNet = teasar( seg; obj_id = one(UInt32) )
    @test length(nodeNet) > 0
end 
