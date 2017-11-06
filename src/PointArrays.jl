module PointArrays

export PointArray

const PointArray = Array{UInt32, 2}

const ZERO_UINT32 = convert(UInt32, 0)
const ONE_UINT32  = convert(UInt32, 1)
const DEFAULT_OFFSET = (ZERO_UINT32, ZERO_UINT32, ZERO_UINT32)

const MAX_BOUNDARY_DISTANCE = 100000

"""
parameter:
    bin_im: binary array. object voxels are false, non-object voxels are true!
"""
function from_binary_image(bin_im::Array{Bool,3})
    x,y,z = findn(!bin_im)
    x = Vector{UInt32}(x)
    y = Vector{UInt32}(y)
    z = Vector{UInt32}(z)
    return hcat(x,y,z)
end 

"""
find points inside an object from a segmentation array. 
"""
function from_seg(seg::Array{T,3}; obj_id::T=convert(T,1),
                     offset::NTuple{3,UInt32}=DEFAULT_OFFSET) where T
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
find out the boundary voxels and represent them as indexes in the point array
"""
function get_boundary_point_indexes(self::Array{T,2}, seg::Array{TSeg,3},
                                    obj_id::TSeg = TSeg(1)) where {T, TSeg}
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
function add_offset!(self::Array{T,2}, offset::NTuple{3,T}) where T
    self[:, 1] .+= offset[1]
    self[:, 2] .+= offset[2]
    self[:, 3] .+= offset[3]
end 

"""
    merge(self::Array{T,3}, other::Array{T,3})
"""
function merge(self::Array{T,2}, other::Array{T,2}) where T
    vcat(self, other)
end 

end # module
