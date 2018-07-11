using BigArrays
using BigArrays.GSDicts
using Base.Test
using RealNeuralNetworks 
using RealNeuralNetworks.Manifests
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs

const cellId = 77497
const VOXEL_SIZE = (80,80,45)

@testset "test basic data structure" begin 
    h = GSDict("gs://neuroglancer/zfish_v1/consensus-20180123/mesh_mip_4"; valueType=Dict{Symbol, Any})
    manifest = h["$(cellId):0"]
    str = split(manifest[:fragments][1], ":")[end]
    @show str
    range = BigArrays.Indexes.string2unit_range( str )
    @test range == [2865:3376, 2017:2528, 16401:16912]
end 

@testset "test manifest iteration" begin
    manifest = Manifest("gs://neuroglancer/zfish_v1/consensus-20180123/mesh_mip_4", 
                        "$(cellId):0", "gs://neuroglancer/zfish_v1/consensus-20180123/40_40_45")
    println("trace a cell...")
    @time nodeNet = Manifests.trace(manifest, cellId)
    swc = SWC( nodeNet )
    SWCs.stretch_coordinates!(swc, VOXEL_SIZE)
    SWCs.save(swc, "/tmp/$(cellId).swc")
end 
