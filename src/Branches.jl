module Branches

include("BoundingBoxes.jl")
# include(joinpath(dirname(@__FILE__), "BoundingBoxes.jl"))
using .BoundingBoxes 

const CLASS = UInt8(0)

export Branch 
type Branch 
    # x,y,z,r
    nodeList    ::Vector{NTuple{4, Float32}}
    class       ::UInt8
    boundingBox ::BoundingBox
end 

function Branch(nodeList::Vector; class=CLASS)
    Branch(nodeList, class, BoundingBox(nodeList))
end 

"""
    Base.length(self::Branch)

the number of nodes contained in this branch 
"""
function Base.length(self::Branch)
    length(self.nodeList) 
end 

"""
split the branch from the node list index to two branches
the indexed node will be included in the first branch 
"""
function Base.split(self::Branch, index::Integer)
    @assert index >=1 && index<length(self)
    nodeList1 = self.nodeList[1:index]     
    nodeList2 = self.nodeList[index+1:end]
    branch1 = Branch(nodeList1; class=self.class)
    branch2 = Branch(nodeList2; class=self.class)
    return branch1, branch2
end 

"""
distance from a point 
"""
function distance_from(self::Branch, point::Tuple)
    distance_from(self, [point[1:3]...])
end 
function distance_from(self::Branch, point::Vector)
    ret = (0,0)
    @assert length(point) == 3 || length(point) == 4
    distance = typemax(Float32)
    for (index, node) in enumerate(nodeList)
        d = norm([node[1:3]...], point[1:3])
        if d < distance
            distance = d
            ret = (d, index)
        end 
    end 
    ret 
end 

function get_bounding_box_distance(self::Branch, point::Union{Tuple, Vector})
    @assert length(point) >= 3
    BoundingBoxes.distance_from(self.boundingBox, point)
end 

end # module
