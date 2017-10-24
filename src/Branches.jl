module Branches

include("BoundingBoxes.jl")
# include(joinpath(dirname(@__FILE__), "BoundingBoxes.jl"))
using .BoundingBoxes 

const CLASS = UInt8(0)
const ZERO_FLOAT32 = Float32(0)

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

###################### properties ###################
function get_node_list(self::Branch) self.nodeList end 
function get_connectivity_matrix( self::Branch ) self.connectivityMatrix end 
function get_bounding_box( self::Branch ) self.boundingBox end 
function get_class( self::Branch ) self.class end 

###################### Base functions ################
function Base.isempty(self::Branch) isempty(self.nodeList) end 

"""
    Base.length(self::Branch)

the number of nodes contained in this branch 
"""
function Base.length(self::Branch)
    length(self.nodeList) 
end 

"""
split the branch from the node list index to two branches
the indexed node will be included in the second branch 
"""
function Base.split(self::Branch, index::Integer)
    @assert index >=1 && index<=length(self)
    nodeList1 = self.nodeList[1:index-1]     
    nodeList2 = self.nodeList[index:end]
    branch1 = Branch(nodeList1; class=self.class)
    branch2 = Branch(nodeList2; class=self.class)
    return branch1, branch2
end

"""
    Base.getindex(self::Branch, index::Integer)
"""
function Base.getindex(self::Branch, index::Integer)
    get_node_list(self)[ index ]
end 

"""
distance from a point 
"""
function distance_from(self::Branch, point::Tuple)
    distance_from(self, [point[1:3]...])
end 
function distance_from(self::Branch, point::Vector)
    @assert !isempty(self)
    ret = (0,0)
    nodeList = get_node_list(self)
    @assert length(point) == 3 || length(point) == 4
    distance = typemax(Float32)
    for (index, node) in enumerate(nodeList)
        d = norm( [node[1:3]...] .- [point[1:3]...] )
        if d < distance
            distance = d
            ret = (d, index)
        end 
    end
    @assert ret[1] > 0
    @assert ret[2] <= length(self)
    @assert ret!=(0,0)
    ret 
end 

function get_bounding_box_distance(self::Branch, point::Union{Tuple, Vector})
    @assert length(point) >= 3
    BoundingBoxes.distance_from(self.boundingBox, point)
end 

function add_offset(self::Branch, offset::Union{Tuple, Vector})
    @assert length(offset) == 3
    nodeList = Vector{NTuple{4,Float32}}()
    for node in self.nodeList
        newNode = map(+, node, [offset..., ZERO_FLOAT32])
        push!(nodeList, newNode)
    end
    Branch(nodeList, self.class, self.boundingBox)    
end 


end # module
