module BoundingBoxes

const ZERO_FLOAT32 = Float32(0)

export BoundingBox 

type BoundingBox
    minCorner   ::NTuple{3,Float32}
    maxCorner   ::NTuple{3,Float32}
end 

function BoundingBox(minCorner::Vector, maxCorner::Vector)
    @assert length(minCorner) == 3
    @assert length(maxCorner) == 3
    BoundingBox((minCorner...), (maxCorner...))
end 

"""
get bounding box from a node list 
"""
function BoundingBox(nodeList::Vector{NTuple{4,Float32}})
    minCorner = [Inf32, Inf32, Inf32]
    maxCorner = [ZERO_FLOAT32, ZERO_FLOAT32, ZERO_FLOAT32]
    for node in nodeList 
        minCorner = min(minCorner, node[1:3])
        maxCorner = max(maxCorner, node[1:3])
    end 
    BoundingBox(minCorner, maxCorner)
end 

function distance_from(self::BoundingBox, point::Union{Tuple, Vector})
    @assert length(point) == 3
    min(norm([self.minCorner...] .- [point...]), 
        norm([self.maxCorner...] .- [point...]))
end

function Base.isequal(self::BoundingBox, other::BoundingBox)
    self.minCorner==other.minCorner && self.maxCorner==other.maxCorner 
end
function Base.:(==)(self::BoundingBox, other::BoundingBox)
    isequal(self, other) 
end 


function Base.union(self::BoundingBox, other::BoundingBox)
    minCorner = map(min, self.minCorner, other.minCorner)
    maxCorner = map(max, self.maxCorner, other.maxCorner)
    BoundingBox(minCorner, maxCorner)
end

end # module
