module PointArrays

const ZERO_UINT32 = convert(UInt32, 0)
const ONE_UINT32  = convert(UInt32, 1)
const DEFAULT_OFFSET = (ZERO_UINT32, ZERO_UINT32, ZERO_UINT32)

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
    return points
end

"""
add offset to points
"""
function add_offset{T}(self::Array{T,3}, offset::NTuple{3,T})
    self[:, 1] += offset[1]
    self[:, 2] += offset[2]
    self[:, 3] += offset[3]
    return self 
end 

"""
    merge(self::Array{T,3}, other::Array{T,3})
"""
function merge{T}(self::Array{T,3}, other::Array{T,3})
    vcat(self, other)
end 

end # module
