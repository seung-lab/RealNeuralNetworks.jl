module PointArrays

const ZERO_UINT32 = convert(UInt32, 0)
const ONE_UINT32  = convert(UInt32, 1)
const DEFAULT_OFFSET = (ZERO_UINT32, ZERO_UINT32, ZERO_UINT32)

const MAX_BOUNDARY_DISTANCE = 100000

"""
find points inside an object from a segmentation array. 
"""
function from_seg{T}(seg::Array{T,3}; obj_id::T=convert(T,1),
                        offset::NTuple{3,UInt32}=DEFAULT_OFFSET)
    # compute the number of object voxels
    nov = 0
    for i in eachindex(seg)
        if seg[i] == obj_id 
            nov += 1
        end 
    end 
    # initialize the points array with offset
    points = Array(UInt32, (nov, 3))
    points[:,1] = offset[1]
    points[:,2] = offset[2]
    points[:,3] = offset[3]
    
    # position of current row
    pos = 0
    for z in ONE_UINT32 : convert(UInt32, size(seg, 3))
        for y in ONE_UINT32 : convert(UInt32, size(seg, 2))
            for x in ONE_UINT32 : convert(UInt32, size(seg, 1))
                if seg[x,y,z] == obj_id
                    pos += 1
                    points[pos,:] .+= [x,y,z]
                end 
            end 
        end 
    end
    # find out the boundary voxel indexes
    boundary_point_indexes = get_boundary_point_indexes(points, seg; obj_id=obj_id)
    return points, boundary_point_indexes 
end

"""
find out the boundary voxels and represent them as indexes in the point array
"""
function get_boundary_point_indexes{T, TSeg}(self::Array{T,2}, seg::Array{TSeg,3};
                                             obj_id::T = T(1))
    # compute the indexes of boundary voxels
    ret = Vector{T}()
    for i in 1:size(self,1)
        point = self[i,:]
        x,y,z = (point...)
        if  z==1 || z==size(seg, 3) ||
            y==1 || y==size(seg, 2) ||
            x==1 || x==size(seg, 1) 
            push!(ret, i)
        else
            if  seg[x-1,y,z]!=obj_id || seg[x+1,y,z]!=obj_id ||
                seg[x,y-1,z]!=obj_id || seg[x,y+1,z]!=obj_id ||
                seg[x,y,z-1]!=obj_id || seg[x,y,z+1]!=obj_id 
                push!(ret, i)
            end 
        end
    end 
    return ret
end 

"""
add offset to points
"""
function add_offset{T}(self::Array{T,2}, offset::NTuple{3,T})
    self[:, 1] += offset[1]
    self[:, 2] += offset[2]
    self[:, 3] += offset[3]
    return self 
end 

"""
    merge(self::Array{T,3}, other::Array{T,3})
"""
function merge{T}(self::Array{T,2}, other::Array{T,2})
    vcat(self, other)
end 

"""
compute Distance from Boundary Field (DBF) based on point cloud and the boundary points 
"""
function compute_DBF{T}( self::Array{T,2}, boundary_point_indexes::Vector )
    num = size(self, 1)
    dbf = Vector{Float32}(num)
    fill!(dbf, Inf32)
    for i in 1:num
        point = self[i,:]
        for bpi in boundary_point_indexes
            boundary = self[bpi, :]
            # filter out some far away boundary points 
            if  abs(point[1]-boundary[1]) < MAX_BOUNDARY_DISTANCE && 
                abs(point[2]-boundary[2]) < MAX_BOUNDARY_DISTANCE && 
                abs(point[3]-boundary[3]) < MAX_BOUNDARY_DISTANCE 
                # compute euclidean distance
                ed = norm(point .- boundary)
                if ed < dbf[i]
                    dbf[i] = ed 
                end 
            end 
        end 
    end
    return dbf
end 

end # module
