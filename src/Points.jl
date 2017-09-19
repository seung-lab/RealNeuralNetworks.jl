module Points
using BigArrays.Chunks

function from_seg{T}(seg::Array{T,3}; obj_id::T=convert(T,1))
    # compute the number of object voxels
    nov = 0
    for i in eachindex(seg)
        if seg[i] == obj_id 
            nov += 1
        end 
    end 
    # initialize the points array
    points = zeros(UInt32, (nov, 3))
    # position of current row
    pos = 0
    for z in UInt32(1):UInt32(size(seg, 3))
        for y in UInt32(1):UInt32(size(seg, 2))
            for x in UInt32(1):UInt32(size(seg, 1))
                if seg[x,y,z] == obj_id
                    pos += 1
                    points[pos,:] = [x,y,z]
                end 
            end 
        end 
    end
    return points
end 

function from_chunk{T}(chk::Chunk; obj_id::T=convert(T,1))
    error("not implemented!")    
end 
end # module
