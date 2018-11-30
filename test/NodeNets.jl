using RealNeuralNetworks  
using Test
using HDF5
using BigArrays
using BigArrays.GSDicts 
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs

using LinearAlgebra

const CELL_ID = UInt32(76880)
const EXPANSION= (UInt32(80), UInt32(80), UInt32(40))
const MIP = 4

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

@testset "test NodeNets module" begin
    seg = create_fake_seg()
    # @time seg = get_seg_from_h5()
    # @time seg = get_seg_from_gs()
    println("building nodeNet ...")
    #@time nodeNet = NodeNet( seg; obj_id = CELL_ID )
    @time nodeNet = NodeNet( seg; obj_id = one(UInt32) )
    NodeNets.add_offset!(nodeNet, (-one(Float32),-one(Float32),-one(Float32)))
    bin = NodeNets.get_neuroglancer_precomputed(nodeNet)
    # open("/tmp/fake.bin", "w") do f write(f, bin)  end 
    @time swc = SWCs.SWC( nodeNet )
    SWCs.stretch_coordinates!(swc, MIP)
    SWCs.save(swc, "/tmp/test.swc")
    @test length(swc) > 1
end 


