module BoundingBoxes

import LinearAlgebra: norm 

const ZERO_FLOAT32 = zero(Float32)

export BoundingBox 

mutable struct BoundingBox
    minCorner   ::NTuple{3,Float32}
    maxCorner   ::NTuple{3,Float32}
end 

function BoundingBox()
    minCorner = (Inf32, Inf32, Inf32)
    maxCorner = (ZERO_FLOAT32, ZERO_FLOAT32, ZERO_FLOAT32)
    BoundingBox( minCorner, maxCorner )
end 

function BoundingBox( minCorner::NTuple{3,Int}, maxCorner::NTuple{3,Int} )
	BoundingBox( 	map(x->convert(Float32,x), minCorner), 
					map(x->convert(Float32,x), maxCorner) )
end 

function BoundingBox( minCorner::Vector{Float32}, maxCorner::Vector{Float32} )
    BoundingBox((minCorner...,), (maxCorner...,))
end 

"""
get bounding box from a node list 
"""
function BoundingBox(nodeList::Union{Vector, Set})
    minCorner = fill(Inf32,3)
    maxCorner = fill(ZERO_FLOAT32, 3)
    for node in nodeList 
        minCorner = map(min, minCorner, node[1:3])
        maxCorner = map(max, maxCorner, node[1:3])
    end 
    BoundingBox(minCorner, maxCorner)
end 

function get_unit_range( self::BoundingBox )
    map((x,y)->floor(Int,x):ceil(Int,y), self.minCorner, self.maxCorner) 
end 

function Base.size(self::BoundingBox)
    map(length, get_unit_range(self))
end 

function Base.isequal(self::BoundingBox, other::BoundingBox)
    self.minCorner==other.minCorner && self.maxCorner==other.maxCorner 
end

function Base.:(==)(self::BoundingBox, other::BoundingBox)
    isequal(self, other) 
end 

"""
    Base.union(self::BoundingBox, other::BoundingBox)
get the bounding box containing these two bounding boxes 
"""
function Base.union(self::BoundingBox, other::BoundingBox)
    minCorner = map(min, self.minCorner, other.minCorner)
    maxCorner = map(max, self.maxCorner, other.maxCorner)
    BoundingBox(minCorner, maxCorner)
end

function isinside(self::BoundingBox, point::Union{Tuple, Vector})
    all( map((x,y)->x>y, self.maxCorner, point[1:3]) ) && 
    all( map((x,y)->x<y, self.minCorner, point[1:3]) )
end 

"""
    to_voxel_space(self::BoundingBox, voxelSize::Union{Tuple, Vector})

"""
function to_voxel_space(self::BoundingBox, voxelSize::Union{Tuple, Vector})
    minCorner = map(fld, self.minCorner, voxelSize)
    maxCorner = map(cld, self.maxCorner, voxelSize)
    BoundingBox(minCorner, maxCorner)
end 

"""
    distance_from(self::BoundingBox, point::Union{Tuple, Vector})

compute the distance from bounding box using a smart way
https://stackoverflow.com/questions/5254838/calculating-distance-between-a-point-and-a-rectangular-box-nearest-point

code in JS
function distance(rect, p) {
  var dx = Math.max(rect.min.x - p.x, 0, p.x - rect.max.x);
  var dy = Math.max(rect.min.y - p.y, 0, p.y - rect.max.y);
  return Math.sqrt(dx*dx + dy*dy);
}
"""
function distance_from(self::BoundingBox, point::Union{Tuple, Vector})
    @assert length(point) == 3 || length(point)==4
    d = map((cmin, p, cmax)-> max(cmin-p, Float32(0), p-cmax), 
                            self.minCorner, point[1:3], self.maxCorner )
    norm([d...])
end


end # module
