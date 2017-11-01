module Branches

include("BoundingBoxes.jl")
# include(joinpath(dirname(@__FILE__), "BoundingBoxes.jl"))
using .BoundingBoxes 

typealias Node NTuple{4,Float64}

const CLASS = UInt8(0)
const ZERO_FLOAT32 = Float32(0)

export Branch 
type Branch 
    # x,y,z,r
    nodeList    ::Vector{Node}
    class       ::UInt8
    boundingBox ::BoundingBox
end 

function Branch(nodeList::Vector; class=CLASS)
    Branch(nodeList, class, BoundingBox(nodeList))
end 

###################### properties ###################
"""
    get_nodes_distance(self::Node, other::Node)
compute the euclidean distance between two nodes 
"""
function get_nodes_distance(self::Node, other::Node)
    norm( [map((x,y)->x-y, self[1:3], other[1:3]) ...])
end 
function get_node_list(self::Branch) self.nodeList end 
function get_connectivity_matrix( self::Branch ) self.connectivityMatrix end 
function get_bounding_box( self::Branch ) self.boundingBox end 
function get_class( self::Branch ) self.class end 

function get_bounding_box_distance(self::Branch, point::Union{Tuple, Vector})
    @assert length(point) >= 3
    BoundingBoxes.distance_from(self.boundingBox, point)
end 

"""
    get_path_length(self::Branch)
accumulate the euclidean distance between neighboring nodes 
"""
function get_path_length(self::Branch)
    ret = 0.0
    for i in 2:length(self)
        ret += norm( [map((x,y)-> x-y, self[i][1:3], self[i-1][1:3] )...] )
    end
    ret
end 

###################### Base functions ################
function Base.endof(self::Branch) length(self) end
function Base.isempty(self::Branch) isempty(self.nodeList) end 

"""
    Base.length(self::Branch)

the number of nodes contained in this branch 
"""
function Base.length(self::Branch)
    length(self.nodeList) 
end 

"""
    Base.merge(self::Branch, other::Branch)
merge two branches  
"""
function Base.merge(self::Branch, other::Branch)
    nodeList1 = get_node_list(self)
    nodeList2 = get_node_list(other)
    mergedNodeList = vcat( nodeList1, nodeList2 )
    # winner taks all!
    class = length(nodeList1)>length(nodeList2)? get_class(self) : get_class(other)
    boundingBox = union( get_bounding_box(self), get_bounding_box(other) )
    Branch(mergedNodeList, class, boundingBox)
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
