using BigArrays
using GSDicts
using Base.Test
using TEASAR 
using TEASAR.Manifests
using TEASAR.NodeNets
using TEASAR.SWCs
using JLD

const cellId = 77497
const VOXEL_SIZE = (80,80,45)

@testset "test basic data structure" begin 
    h = GSDict("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4"; valueType=Dict{Symbol, Any})
    manifest = h["770048087:0"]
    str = split(manifest[:fragments][1], ":")[end]
    range = BigArrays.Indexes.string2unit_range( str )
    @test range == [2969:3480, 1777:2288, 16913:17424]
end 

@testset "test manifest iteration" begin
    manifest = Manifest("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4", 
                        "$(cellId):0", "gs://neuroglancer/zfish_v1/consensus-20170829/80_80_45")
    nodeNet = Manifests.trace(manifest, cellId)
    swc = SWC( nodeNet )
    SWCs.stretch_coordinates!(swc, VOXEL_SIZE)
    SWCs.save(swc, "/tmp/$(cellId).swc")
end 
