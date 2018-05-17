module Segments

using RealNeuralNetworks.Utils.BoundingBoxes
include("Synapses.jl")
using .Synapses

const Node = NTuple{4,Float32}
const SynapseList = SparseVector{Synapse, Int}

const CLASS = zero(UInt8)

export Segment 
mutable struct Segment 
    # list of tuple (x,y,z,r)
    nodeList        ::Vector{Node}
    class           ::UInt8
    boundingBox     ::BoundingBox
    preSynapseList  ::SynapseList 
    postSynapseList ::SynapseList 
end 

function Segment(nodeList::Vector{Node}; 
                 class::UInt8=CLASS, boundingBox=BoundingBox(nodeList),
                 preSynapseList::SynapseList  = spzeros(Synapse, length(nodeList)),
                 postSynapseList::SynapseList = spzeros(Synapse, length(nodeList)))
    Segment(nodeList, class, boundingBox, preSynapseList, postSynapseList)
end 

###################### properties ###################
"""
    get_nodes_distance(self::Node, other::Node)
compute the euclidean distance between two nodes 
"""
@inline function get_nodes_distance(self::Union{Vector,Tuple}, other::Union{Vector,Tuple})
    norm( [map((x,y)->x-y, self[1:3], other[1:3]) ...])
end 
@inline function get_node_list(self::Segment) self.nodeList end 
@inline function get_connectivity_matrix( self::Segment ) self.connectivityMatrix end 
@inline function get_bounding_box( self::Segment ) self.boundingBox end 
@inline function get_class( self::Segment ) self.class end 
@inline function get_pre_synapse_list( self::Segment ) self.preSynapseList end 
@inline function get_post_synapse_list( self::Segment ) self.postSynapseList end
@inline function get_pre_synapse( self::Segment, index::Int ) self.preSynapseList[index] end
@inline function get_post_synapse( self::Segment, index::Int ) 
                                        self.postSynapseList[index] end

@inline function get_bounding_box_distance(self::Segment, point::Union{Tuple, Vector})
    @assert length(point) >= 3
    BoundingBoxes.distance_from(self.boundingBox, point)
end 

"""
    get_path_length(self::Segment; nodeIndex::Int=length(self))
accumulate the euclidean distance between neighboring nodes 
"""
@inline function get_path_length(self::Segment; nodeIndex::Int=length(self))
    ret = 0.0
    for i in 2:nodeIndex
        ret += norm( [map((x,y)-> x-y, self[i][1:3], self[i-1][1:3] )...] )
    end
    ret
end

@inline function get_radius_list( self::Segment ) map(n->n[4], self) end 

"""
    get_tail_head_radius_ratio( self::Segment )
the spine is normally thick in tail, and thin in the head. 
ratio = max_tail / mean_head
The head should point to dendrite. This is a very good feature to identify spine.
"""
@inline function get_tail_head_radius_ratio( self::Segment )
    radiusList = get_radius_list( self )
    N = length(self)
    headRadiusList = radiusList[1:cld(N,2)]
    tailRadiusList = radiusList[cld(N,2):N]
    maximum(tailRadiusList) / mean(headRadiusList)
end 

"""
    get_tortuosity( self::Segment )
the ratio of the actual path length to the euclidean distance between head and tail node 
"""
@inline function get_tortuosity(self::Segment)
    if length(self) == 1 
        return 1.0
    end 
    pathLength = get_path_length(self)
    euclideanLength = get_nodes_distance( self[1], self[end] )
    @assert self[1]!=self[end] "segment start is the same with the end: $(self)"
    @assert euclideanLength != 0.0
    pathLength / euclideanLength 
end 

@inline function get_center(nodeList::Vector{Node})
    (map(i->mean(y->y[i], nodeList), 1:4)...)
end

@inline function get_center(self::Segment, range::UnitRange)
    center = get_center( self[range] )
end 

###################### Base functions ################
function Base.start( self::Segment ) 1 end 
function Base.next( self::Segment, state::Integer ) get_node_list(self)[state], state+1 end 
function Base.done( self::Segment, state::Integer ) state > length(self) end 

function Base.endof(self::Segment) length(self) end
function Base.isempty(self::Segment) isempty(self.nodeList) end 

"""
    Base.length(self::Segment)

the number of nodes contained in this segment 
"""
function Base.length(self::Segment)
    length(self.nodeList) 
end 

"""
    Base.merge(self::Segment, other::Segment)
merge two segmentes  
"""
function Base.merge(self::Segment, other::Segment)
    nodeList1 = get_node_list(self)
    nodeList2 = get_node_list(other)
    mergedNodeList = vcat( nodeList1, nodeList2 )
    # winner taks all!
    class = length(nodeList1)>length(nodeList2) ? get_class(self) : get_class(other)
    boundingBox = union( get_bounding_box(self), get_bounding_box(other) )
    Segment(mergedNodeList; class=class, boundingBox=boundingBox)
end 

"""
split the segment from the node list index to two segmentes
the indexed node will be included in the second segment 
"""
function Base.split(self::Segment, index::Integer)
    @assert index >=1 && index<=length(self)
    local nodeList1::Vector{NTuple{4, Float32}}
    local nodeList2::Vector{NTuple{4, Float32}}
    if index==1
        if length(self)==1
            nodeList1 = self.nodeList[1]
            nodeList2 = self.nodeList[1]
        else 
            nodeList1 = self.nodeList[1:index]     
            nodeList2 = self.nodeList[index+1:end]
        end 
    else 
        nodeList1 = self.nodeList[1:index-1]     
        nodeList2 = self.nodeList[index:end]
    end 
    segment1 = Segment(nodeList1; class=self.class)
    segment2 = Segment(nodeList2; class=self.class)
    @assert length(nodeList1) > 0
    @assert length(nodeList2) > 0
    return segment1, segment2
end

"""
    Base.getindex(self::Segment, index::Integer)
"""
@inline function Base.getindex(self::Segment, index::Integer)
    get_node_list(self)[ index ]
end

@inline function Base.getindex(self::Segment, range::UnitRange)
    get_node_list(self)[range]
end 

"""
distance from a point 
"""
function distance_from(self::Segment, point::Tuple)
    distance_from(self, [point[1:3]...])
end 
function distance_from(self::Segment, point::Vector)
    @assert !isempty(self)
    ret = (zero(Float32), zero(Int))
    nodeList = get_node_list(self)
    @assert !isempty(nodeList)
    @assert length(point) == 3 || length(point) == 4
    distance = typemax(Float32)
    for (index, node) in enumerate(nodeList)
        d = norm( [node[1:3]...] .- [point[1:3]...] )
        if d < distance
            distance = d
            ret = (d, index)
        end 
    end
    @assert ret[2] <= length(self)
    ret 
end 

################## manipulation ###############################

function attach_pre_synapse!(self::Segment, nodeIndex::Int, synapse::Synapse)
    self.preSynapseList[ nodeIndex ] = synapse 
end 

function attach_post_synapse!(self::Segment, nodeIndex::Int, synapse::Synapse)
    self.postSynapseList[ nodeIndex ] = synapse
end 

function add_offset(self::Segment, offset::Union{Tuple, Vector})
    @assert length(offset) == 3
    nodeList = Vector{NTuple{4,Float32}}()
    for node in self.nodeList
        newNode = map(+, node, [offset..., zero(Float32)])
        push!(nodeList, newNode)
    end
    Segment(nodeList, self.class, self.boundingBox)    
end 

function remove_node(self::Segment, removeNodeIndex::Integer)
    newNodeList = Vector{NTuple{4, Float32}}()
    for (index,node) in enumerate(get_node_list(self))
        if index != removeNodeIndex 
            push!(newNodeList, node)
        end 
    end 
    Segment(newNodeList; class = get_class(self))
end 

function remove_redundent_nodes!(self::Segment)
    nodeList = get_node_list(self)
    newNodeList = Vector{NTuple{4, Float32}}()
    for index in 1:length(nodeList)-1
        if nodeList[index] != nodeList[index+1]
            push!(newNodeList, nodeList[index])
        end 
    end 
    # get the final node
    push!(newNodeList, nodeList[end])
    self.nodeList = newNodeList 
end 

end # module
