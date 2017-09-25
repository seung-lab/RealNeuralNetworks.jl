using BigArrays
using GSDicts
using Base.Test

@testset "test manifest" begin 
    h = GSDict("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4"; valueType=Dict{Symbol, Any})
    manifest = h["770048087:0"]
    str = split(manifest[:fragments][1], ":")[end]
    range = BigArrays.Indexes.string2unit_range( str )
    @test range == [2968:3480, 1776:2288, 16912:17424]
end 
