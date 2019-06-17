using BigArrays
using BigArrays.BinDicts
using Test
using RealNeuralNetworks 
using RealNeuralNetworks.Manifests
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs
using JSON 

const cellId = 77497
const VOXEL_SIZE = (80,80,45)

@testset "test basic data structure" begin 
    h = BinDict("/neuroglancer/zfish_v1/consensus-20180123/mesh_mip_4")
    manifest = String(h["$(cellId):0"])
    manifest = JSON.parse(manifest; dicttype=Dict{Symbol, Any})
    
    str = split(manifest[:fragments][1], ":")[end]
    @show str
    range = BigArrays.Indexes.string2unit_range( str )
    @test range == [2865:3376, 2017:2528, 16401:16912]
end 

@testset "test manifest iteration" begin
    manifest = Manifest("/neuroglancer/zfish_v1/consensus-20180123/mesh_mip_4", 
                        "$(cellId):0", "/neuroglancer/zfish_v1/consensus-20180123/40_40_45")
    println("trace a cell...")
    @time nodeNet = Manifests.trace(manifest, cellId)
    swc = SWC( nodeNet )
    SWCs.stretch_coordinates!(swc, VOXEL_SIZE)
    SWCs.save(swc, "/tmp/$(cellId).swc")
end 
