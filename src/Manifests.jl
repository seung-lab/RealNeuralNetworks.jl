module Manifests
include("DBFs.jl"); using .DBFs;
include("PointArrays.jl"); using .PointArrays;

using RealNeuralNetworks

# using JLD2
using BigArrays
using BigArrays.GSDicts

using OffsetArrays
using JSON
#import Distributed: pmap 

const MIP_LEVEL = 4

export Manifest

struct Manifest
    # the bigarray for cutout of segmentation 
    ba          ::AbstractBigArray
    # the id of object
    obj_id      ::Integer
    # the unit range list for cutout
    rangeList   ::Vector
end 

"""
Parameters:
    dir 
"""
function Manifest(manifestDirPath::AbstractString, manifestKey::AbstractString, 
                  bigArrayPath::AbstractString, mip::Integer )
    #println("big array path: ", bigArrayPath)
    #println("manifestDirPath: ", manifestDirPath)

    ba = BigArray( GSDict( bigArrayPath ); mip=mip, mode=:multithreads )
    h = GSDict( manifestDirPath )
    ranges_str_list = h[manifestKey][:fragments]
    # @show ranges_str_list

    obj_id = Meta.parse( split(ranges_str_list[1], ":")[1] )
    obj_id = convert(eltype(ba), obj_id)
    ranges_str_list = map(x-> split(x,":")[end], ranges_str_list)
    rangesList = map( BigArrays.Indexes.string2unit_range, ranges_str_list )
    # convert from z,y,x to x,y,z
    # this is due to a bug in chunkflow, we should not need this in the future
    rangesList = map( reverse, rangesList)
    # @show rangesList 
    Manifest( ba, obj_id, rangesList )
end

function Base.length(self::Manifest) length(get_range_list(self)) end 

function get_range_list(self::Manifest) self.rangeList end 

"""
the voxel offset in the neuroglancer precomputed info file 
"""
function get_voxel_offset(self::Manifest)
    voxel_offset =  self.ba.kvStore.configDict[:offset]
    #real voxel offset should based on the highest resolution 5x5x45
    # we are using mip level 4 with (80x80x45 nm) now 
    # voxel_offset = map((x,y)->UInt32(x*y), voxel_offset, (2^MIP_LEVEL, 2^MIP_LEVEL,1))
    voxel_offset = Vector{UInt32}(voxel_offset)
    @show voxel_offset
    return (voxel_offset...,)
end 

"""
iterate the chunks containing the neuron with specified neuronId
build point cloud and dbf when iterating the chunks 
"""
function trace(self::Manifest)
    println("extract point clouds and distance from boundary fields ...")
    
    pointCloudDBFList = map(x->_get_point_cloud_dbf(self, x), self.rangeList )

    pointClouds = map( x->x[1], pointCloudDBFList )
    pointCloud = vcat(pointClouds ...)
    dbfs = map(x->x[2], pointCloudDBFList)
    dbf = vcat(dbfs ...)
    # save temporal variables for debug
    # @save "/tmp/$(neuronId).jld" pointClouds, pointCloud, dbf
    println("skeletonization from global point cloud and dbf using RealNeuralNetworks algorithm...")
    @time neuron = teasar(pointCloud; dbf=dbf) 
    return neuron
end 

function _get_point_cloud_dbf(self::Manifest, ranges::Vector)
    # example: [2456:2968, 1776:2288, 16400:16912]
    offset = (map(x-> UInt32(x.start-1), ranges)...,)
    seg = self.ba[ranges...] |> parent
    bin_im = DBFs.create_binary_image( seg; obj_id = self.obj_id )
    @assert !all(bin_im)
    point_cloud = PointArrays.from_binary_image( bin_im )
    @assert !isempty(point_cloud)
    # distance from boundary field
    dbf = DBFs.compute_DBF(point_cloud, bin_im)
    PointArrays.add_offset!(point_cloud, offset)
    # no need to use voxel_offset since the file name encoded the global coordinate
    # PointArrays.add_offset!(point_cloud, get_voxel_offset(self)) 
    return point_cloud, dbf
end  

end # module
