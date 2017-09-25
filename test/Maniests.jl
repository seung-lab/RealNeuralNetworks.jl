using BigArrays
using GSDicts
using Base.Test
using TEASAR 
using TEASAR.Manifests
using JLD

@testset "test basic data structure" begin 
    h = GSDict("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4"; valueType=Dict{Symbol, Any})
    manifest = h["770048087:0"]
    str = split(manifest[:fragments][1], ":")[end]
    range = BigArrays.Indexes.string2unit_range( str )
    @test range == [2969:3480, 1777:2288, 16913:17424]
end 

@testset "test manifest iteration" begin
    manifest = Manifest("gs://neuroglancer/zfish_v1/consensus-20170829/mesh_mip_4", 
                        "76880:0", "gs://neuroglancer/zfish_v1/consensus-20170829/80_80_45")
    pointCloudDBFList = pmap( identity, manifest )
    pointClouds = map( x->x[1], pointCloudDBFList )
    pointCloud = vcat(pointClouds ...)
    dbfs = map(x->x[2], pointCloudDBFList)
    dbf = vcat(dbfs ...)
    save("/tmp/point_clouds.jld", "point_clouds", pointClouds, 
            "point_cloud", pointCloud, "dbf", dbf)
    swc = skeletonize(pointCloud; dbf=dbf) 
    TEASAR.SWCs.save(swc, "/tmp/76880.swc")
end 
