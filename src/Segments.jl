module Segments

import LinearAlgebra: norm, dot
import Statistics: mean, std 

using RealNeuralNetworks.Utils.BoundingBoxes
include("Synapses.jl"); using .Synapses

const Node = NTuple{4,Float32}

# TO-DO: make the data type more general
const SynapseList = Vector{Union{Missing, Vector{Synapse{Float32}}}}

# classes following SWC format 
const AXON_CLASS = UInt8(2)
const DENDRITE_CLASS = UInt8(3)
const SOMA_CLASS = one(UInt8)
const UNDEFINED_CLASS = zero(UInt8)

export Segment 
mutable struct Segment{T}  
    # list of tuple (x,y,z,r)
    nodeList        ::Vector{NTuple{4,T}}
    class           ::UInt8
    #boundingBox     ::BoundingBox
    preSynapseList  ::SynapseList 
    postSynapseList ::SynapseList 
end 

function Segment()
    nodeList = Vector{Node}()
    Segment(nodeList)
end

function Segment(nodeList::Vector{Node}; 
                 class::UInt8=UNDEFINED_CLASS,
                 preSynapseList::Vector{X}  = SynapseList(missing, length(nodeList)),
                 postSynapseList::Vector{X} = SynapseList(missing, length(nodeList))) where {X<:Union{Missing, Vector{Synapse{Float32}}}}
    Segment(nodeList, class, preSynapseList, postSynapseList)
end 

###################### IO ###########################


###################### properties ###################

"""
    get_features(self::Segment) 
aggregate features to a named tuple 
"""
function get_features(self::Segment)
    (pathLength=get_path_length(self), 
        surfaceArea=get_surface_area(self), 
        volume=get_volume(self),
        meanRadius=mean(get_radius_list(self)),
        stdRadius=std(get_radius_list(self)),
        numPreSynapses=get_num_pre_synapses(self),
        numPostSynapses=get_num_post_synapses(self),
        tortuosity=get_tortuosity(self))
end  

"""
    get_nodes_distance(self::Node, other::Node)
compute the euclidean distance between two nodes 
"""
@inline function get_nodes_distance(self::Union{Vector,Tuple}, other::Union{Vector,Tuple})
    norm( [map((x,y)->x-y, self[1:3], other[1:3]) ...])
end 

@inline function get_node_list(self::Segment) self.nodeList end 
@inline function get_connectivity_matrix( self::Segment ) self.connectivityMatrix end 
@inline function get_bounding_box( self::Segment ) BoundingBox(get_node_list(self)) end 
@inline function get_class( self::Segment ) self.class end 
@inline function get_pre_synapse_list( self::Segment ) self.preSynapseList end 
@inline function get_post_synapse_list( self::Segment ) self.postSynapseList end
@inline function get_pre_synapse( self::Segment, index::Int ) self.preSynapseList[index] end
@inline function get_post_synapse( self::Segment, index::Int ) self.postSynapseList[index] end

@inline function get_bounding_box_distance(self::Segment, point::Union{Tuple, Vector})
    @assert length(point) >= 3
    boundingBox = get_bounding_box(self) 
    BoundingBoxes.distance_from(boundingBox, point)
end 


@inline function euclidean_distance( n1::NTuple{3,T}, n2::NTuple{3,T}) where T
    norm([map((x,y)->x-y, n1, n2)...])  
end 

"""
    get_path_length(self::Segment; nodeId::Int=length(self))
accumulate the euclidean distance between neighboring nodes 
"""
@inline function get_path_length(self::Segment; nodeId::Int=length(self))
    ret = 0.0
    for i in 2:nodeId
        ret += euclidean_distance(self[i][1:3], self[i-1][1:3])
    end
    ret
end

@inline function get_num_pre_synapses(self::Segment)
    synapseList = get_pre_synapse_list(self)
    n = 0
    for synapses in synapseList
        if !ismissing(synapses)
            # there exists some synapses
            n += length(synapses)
        end
    end
    return n
end 

@inline function get_num_post_synapses(self::Segment)
    synapseList = get_post_synapse_list(self)
    n = 0
    for synapses in synapseList
        if !ismissing(synapses)
            # there exists some synapses
            n += length(synapses)
        end
    end
    return n
end 

"""
    get_pre_synapse_density(self::Segment)
note that the unit is # / micron 
"""
@inline function get_pre_synapse_density(self::Segment)
    numPreSynapses = get_num_pre_synapses(self)
    pathLength = get_path_length(self)
    numPreSynapses / pathLength * 1000
end 

"""
    get_post_synapse_density(self::Segment)
note that the unit is # / micron 
"""
@inline function get_post_synapse_density(self::Segment)
    numPostSynapses = get_num_post_synapses(self)
    pathLength = get_path_length(self)
    numPostSynapses / pathLength * 1000
end 


@inline function get_radius_list( self::Segment ) map(n->n[4], self) end 

"""
    get_surface_area(self::Segment)
frustum-based:  
http://www.analyzemath.com/Geometry_calculators/surface_volume_frustum.html
"""
function get_surface_area(self::Segment{T}) where T 
    ret = zero(T)  
    for i in 2:length(self) 
        # average diameter
        r1 = self[i][4]
        r2 = self[i-1][4]
        h = euclidean_distance(self[i][1:3], self[i-1][1:3])
        ret += pi* (r1+r2) * sqrt(h*h + (r1-r2)*(r1-r2))
    end
    ret::T 
end

"""
    get_volume(self::Segment)
compute frustum-based volume 
http://jwilson.coe.uga.edu/emt725/Frustum/Frustum.cone.html 
"""
function get_volume(self::Segment{T}) where T
    ret = zero(T) 
    for i in 2:length(self)
        r1 = self[i-1][4]
        r2 = self[i][4]
        h = euclidean_distance(self[i-1][1:3], self[i][1:3])
        ret += pi * h * (r1*r1 + r1*r2 + r2*r2) / T(3)
    end 
    ret 
end 

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
    (map(i->mean(y->y[i], nodeList), 1:4)...,)
end

@inline function get_center(self::Segment, range::UnitRange)
    center = get_center( self[range] )
end 

###################### Base functions ################

function Base.show(self::Segment)
    if get_class(self) == AXON_CLASS 
        class = "axon"
    elseif get_class(self) == DENDRITE_CLASS 
        class = "dendrite"
    else 
        class = "unknown"
    end   
    println("segment is a ", class, " and have ", length(self), " nodes, ", 
            get_num_pre_synapses(self), " presynapses, ", 
            get_num_post_synapses(self), " postsynapses.")
    nothing 
end 

function Base.iterate(self::Segment, state::Int=1)
    if state > length(self) 
        return nothing 
    else 
        return get_node_list(self)[state], state+1  
    end  
end 

@inline function Base.lastindex(self::Segment) length(self) end 
@inline function Base.isempty(self::Segment) isempty(self.nodeList) end 

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
    # merge synapse list
    preSynapseList1 = get_pre_synapse_list(self)
    preSynapseList2 = get_pre_synapse_list(other)
    postSynapseList1 = get_post_synapse_list(self)
    postSynapseList2 = get_post_synapse_list(other)

    mergedPreSynapseList = vcat( preSynapseList1, preSynapseList2 )
    mergedPostSynapseList = vcat( postSynapseList1, postSynapseList2 )

    Segment(mergedNodeList; class=class, preSynapseList=mergedPreSynapseList, 
            postSynapseList=mergedPostSynapseList)
end 

"""
split the segment from the node list index to two segmentes
the indexed node will be included in the second segment 
"""
function Base.split(self::Segment{T}, index::Integer) where T
    @assert index >=1 && index<=length(self)
    idx1 = index - 1
    idx2 = index 

    if index==1
        if length(self)==1
            idx1 = 1
            idx2 = 1
        else
            idx1 = index
            idx2 = index + 1 
        end 
    end
    nodeList1 = self.nodeList[1:idx1]
    nodeList2 = self.nodeList[idx2:end]
    preSynapseList1 = self.preSynapseList[1:idx1]
    preSynapseList2 = self.preSynapseList[idx2:end]
    postSynapseList1 = self.postSynapseList[1:idx1]
    postSynapseList2 = self.postSynapseList[idx2:end]

    segment1 = Segment(nodeList1; class=self.class, preSynapseList=preSynapseList1, postSynapseList=postSynapseList1)
    segment2 = Segment(nodeList2; class=self.class, preSynapseList=preSynapseList2, postSynapseList=postSynapseList2)
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
function distance_from(self::Segment{T}, point::Vector) where T
    @assert !isempty(self)
    ret = (zero(T), zero(Int))
    nodeList = get_node_list(self)
    @assert !isempty(nodeList)
    @assert length(point) == 3 || length(point) == 4
    distance = typemax(T)
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

"""
    attach_pre_synapse!(self::Segment, nodeId::Int, synapse::Synapse{T})
"""
@inline function attach_pre_synapse!(self::Segment, nodeId::Int, synapse::Synapse{T}) where T
    if ismissing(self.preSynapseList[ nodeId ])
        self.preSynapseList[ nodeId ] = [synapse]
    else
        push!(self.preSynapseList[nodeId], synapse)
    end
    nothing 
end 

"""
    attach_post_synapse!(self::Segment, nodeId::Int, synapse::Synapse)
"""
@inline function attach_post_synapse!(self::Segment, nodeId::Int, synapse::Synapse{T}) where T
    if ismissing(self.postSynapseList[nodeId])
        self.postSynapseList[ nodeId ] = [synapse]
    else
        push!(self.postSynapseList[nodeId], synapse)
    end
    nothing
end 

function adjust_class!(self::Segment) 
    numPreSynapses = get_num_pre_synapses(self)
    numPostSynapses = get_num_post_synapses(self)
    if self.class == UNDEFINED_CLASS 
        if numPreSynapses > numPostSynapses
            # mostly presynapses, so this is an axon 
            # Note that this simple metric fails in the axonal hillock, 
            # where there exist a lot of post synapses.
            self.class = AXON_CLASS 
        elseif numPreSynapses < numPostSynapses
            # mostly postsynapses, so this is a dendrite 
            self.class = DENDRITE_CLASS   
        elseif numPostSynapses == 0 && get_path_length(self)>5000
            # a long segment without both pre and post synapses was considered axon 
            self.class = AXON_CLASS
        else 
            self.class = DENDRITE_CLASS 
        end 
    end 
end 

function add_offset(self::Segment{T}, offset::Union{Tuple, Vector}) where T
    @assert length(offset) == 3
    nodeList = Vector{NTuple{4,T}}()
    for node in self.nodeList
        newNode = map(+, node, [offset..., zero(T)])
        push!(nodeList, newNode)
    end
    Segment(nodeList, self.class)    
end 

@inline function remove_node(self::Segment, removeId::Int)
    remove_nodes(self, removeId:removeId)
end 

"""
    remove_nodes(self::Segment, removeIdRange::UnitRange{Int})
remove nodes from a segment
"""
function remove_nodes(self::Segment{T}, removeIdRange::UnitRange{Int}) where T  
    newLength = length(self) - length(removeIdRange)
    @assert newLength >= 0
    if newLength == 0 
        return Segment()
    end 
    newNodeList = Vector{NTuple{4, T}}()
    sizehint!(newNodeList, newLength)

    for (index,node) in enumerate(get_node_list(self))
        if !(index in removeIdRange)
            push!(newNodeList, node)
        end 
    end
    
    preSynapseList = get_pre_synapse_list(self)
    postSynapseList = get_post_synapse_list(self)
    newPreSynapseList = SynapseList(missing, newLength) 
    newPostSynapseList = SynapseList(missing, newLength)

    newPreSynapseList[1:removeIdRange.start-1]  = preSynapseList[1:removeIdRange.start-1]
    newPostSynapseList[1:removeIdRange.start-1] = postSynapseList[1:removeIdRange.start-1]
    
    newPreSynapseList[removeIdRange.start:end]  = preSynapseList[removeIdRange.stop+1:end]
    newPostSynapseList[removeIdRange.start:end] = postSynapseList[removeIdRange.stop+1:end]

    Segment(newNodeList; class = get_class(self), 
            preSynapseList=newPreSynapseList, postSynapseList=newPostSynapseList)
end 

"""
    remove_redundent_nodes!(self::Segment)
remove neighboring nodes that is the same. 
"""
function remove_redundent_nodes!(self::Segment{T}) where T
    nodeList = get_node_list(self)
    newNodeList = Vector{NTuple{4, T}}()
    newPreSynapseList = SynapseList()
    newPostSynapseList = SynapseList()
    for index in 1:length(nodeList)-1
        if nodeList[index] != nodeList[index+1]
            push!(newNodeList, nodeList[index])
            push!(newPreSynapseList, self.preSynapseList[index])
            push!(newPostSynapseList, self.postSynapseList[index])
        else
            # merge synapses
            newPreSynapseList[index] = (newPreSynapseList[index]..., self.preSynapseList[index+1]...)
            newPostSynapseList[index] = (newPostSynapseList[index]..., self.postSynapseList[index+1]...)
        end 
    end 
    # get the final node
    push!(newNodeList, nodeList[end])
    self.nodeList = newNodeList 
end 

end # module
