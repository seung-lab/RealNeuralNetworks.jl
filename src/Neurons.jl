module Neurons
include("Segments.jl")
using .Segments
using ..NodeNets
using ..SWCs
using ..Utils.BoundingBoxes 
using .Segments.Synapses 

using LsqFit
using ImageFiltering
using OffsetArrays
using DataFrames
using SparseArrays
#import LinearAlgebra: norm, dot, normalize 
using LinearAlgebra 
using Serialization  
using Statistics 
using NearestNeighbors

const EXPANSION = (one(UInt32), one(UInt32), one(UInt32))
const VOXEL_SIZE = (500, 500, 500)

export Neuron

mutable struct Neuron{T} 
    # x,y,z,r, the coordinates should be in physical coordinates
    segmentList ::Vector{Segment{T}}
    connectivityMatrix ::SparseMatrixCSC{Bool, Int}
end 

"""
    Neuron
a neuron modeled by interconnected segments 
"""
function Neuron(nodeNet::NodeNet{T}) where T
    # the properties from nodeNet
    nodes = NodeNets.get_node_list(nodeNet)
    radii = NodeNets.get_radii(nodeNet)
    # TO-DO: add node class information
    #nodeClassList = NodeNets.get_node_class_list(nodeNet)
    # flags labeling whether this node was collected to the net
    collectedFlagVec = falses(length(nodes))
    # connectivity matrix of nodes in nodeNet 
    nodesConnectivityMatrix  = NodeNets.get_connectivity_matrix(nodeNet)
    # locate the root node with largest radius
    # theoritically this should be the center of soma 
    _, rootNodeId = findmax(radii)
    seedNodeIdList::Vector = [rootNodeId]
    # grow the main net
    neuron = Neuron!(rootNodeId, nodeNet, collectedFlagVec)
    # println(sum(collectedFlagVec), " visited voxels in ", length(collectedFlagVec))
    #return neuron # only use the main segment for evaluation

    while !all(collectedFlagVec)
        # println(sum(collectedFlagVec), " visited voxels in ", length(collectedFlagVec))
        # there exist some uncollected nodes 
        # find the uncollected node that is closest to the terminal nodes as seed
        # terminal node should have high priority than normal nodes. 
        seedNodeId2 = find_seed_node_id(neuron, nodeNet, collectedFlagVec)

        mergingSegmentId1, mergingNodeIdInSegment1, weightedTerminalDistance = 
                                find_merging_terminal_node_id(neuron, nodeNet[seedNodeId2])
        # for some spine like structures, they should connect to nearest node
        # rather than terminal point. 
        closestSegmentId1, closestNodeIdInSegment1, closestDistance = 
                                find_closest_node(neuron, nodeNet[seedNodeId2][1:3])
        if closestDistance < weightedTerminalDistance / T(2)
            mergingSegmentId1 = closestSegmentId1 
            mergingNodeIdInSegment1 = closestNodeIdInSegment1 
        end 

        # grow a new net from seed, and mark the collected nodes 
        subnet = Neuron!(seedNodeId2, nodeNet, collectedFlagVec) 
        # merge the subnet to main net
        #@show get_num_nodes(neuron)
        #@show get_num_nodes(subnet)
        totalNumNodes = get_num_nodes(neuron) + get_num_nodes(subnet)
        neuron = merge( neuron, subnet, 
                        mergingSegmentId1, mergingNodeIdInSegment1)
        @assert totalNumNodes == get_num_nodes(neuron)
        #@show get_num_nodes(neuron)
    end
    @assert length(nodeNet) == get_num_nodes(neuron) 
                    "$(length(nodeNet)) !== $(get_num_nodes(neuron))"
    neuron     
end 

"""
    Neuron!(seedNodeId::Integer, nodeNet::NodeNet, collectedFlagVec::BitArray{1})

build a net from a seed node using connected component
mark the nodes in this new net as collected, so the collectedFlagVec was changed.
"""
function Neuron!(seedNodeId::Integer, nodeNet::NodeNet{T}, collectedFlagVec::BitArray{1}) where T
    # initialization
    segmentList = Vector{Segment{T}}()

    parentSegmentIdList = Vector{Int}()
    childSegmentIdList  = Vector{Int}()

    nodes = NodeNets.get_node_list(nodeNet)
    nodesConnectivityMatrix = NodeNets.get_connectivity_matrix( nodeNet )
    @assert length(collectedFlagVec) == length(nodes)
    
    # the seed of segment should record both seed node index 
    # and the the parent segment index
    # the parent segment index of root segment is -1 
    segmentSeedList = [(seedNodeId, -1)]
    # depth first search
    while !isempty(segmentSeedList)
        seedNodeId, segmentParentId = pop!(segmentSeedList)
        # grow a segment from seed
        nodeListInSegment = Vector{NTuple{4, T}}()
        seedNodeIdList = [seedNodeId]
        while true
            # construct this segment
            seedNodeId = pop!(seedNodeIdList)
            # push the seed node
            push!(nodeListInSegment, nodes[seedNodeId])
            # label this node index as collected
            collectedFlagVec[ seedNodeId ] = true 
            # println(sum(collectedFlagVec), " visited voxels in ", length(collectedFlagVec))

            # find the connected nodes
            connectedNodeIdList,_ = findnz(nodesConnectivityMatrix[:, seedNodeId])
            # exclude the collected nodes
            connectedNodeIdList = connectedNodeIdList[ .!collectedFlagVec[connectedNodeIdList] ]

            if length(connectedNodeIdList) == 1
                # belong to the same segment
                push!(seedNodeIdList, connectedNodeIdList[1])
            else
                # finish constructing this segment
                # because this is the terminal branching point or multiple branching points
                segment = Segment(nodeListInSegment)
                push!(segmentList, segment)
                if segmentParentId != -1
                    # this is not the root segment, establish the segment connection
                    # @assert length(segment) > 1
                    push!(parentSegmentIdList, segmentParentId)
                    push!(childSegmentIdList,  length(segmentList))
                end 
                # seed new segments
                # if this is terminal segment, no seed will be pushed
                for connectedNodeId in connectedNodeIdList
                    # the length of segmentList will be the parent segment ID of the seeded segment
                    parentSegmentId = length(segmentList)
                    push!(segmentSeedList, (connectedNodeId, parentSegmentId))
                end 
                break
            end 
        end 
    end
    
    for segment in segmentList 
        @assert !isempty(segment)
    end
    @assert length(parentSegmentIdList) == length(childSegmentIdList)

    # note that the connectivity matrix should be a square matrix for easier use
    connectivityMatrix = sparse(parentSegmentIdList, childSegmentIdList, true, 
                                length(segmentList), length(segmentList))
    
    for childSegmentId in 1:length(childSegmentIdList)
        parentSegmentIdList,_ = findnz(connectivityMatrix[:, childSegmentId])
        # segment should only have one parent 
        @assert length(parentSegmentIdList) <= 1 
    end 
    Neuron{T}(segmentList, connectivityMatrix)
end 

function Neuron( seg::Array{ST,3}; obj_id::ST = convert(ST,1), 
                     expansion::NTuple{3,UInt32}=EXPANSION ) where ST
    nodeNet = NodeNet( seg; obj_id = obj_id, expansion = expansion )
    Neuron( nodeNet )
end

function Neuron( swc::SWC )
    nodeNet = NodeNet( swc )
    # respect the point class 
    Neuron( nodeNet )
end 

"""
    Neuron( swcString::AbstractString )
transform string from swc content to Neuron 
"""
function Neuron( swcString::AbstractString )
    swc = SWC(swcString)
    Neuron( swc )
end

######################### IO ################
function load(fileName::AbstractString)
    if endswith(fileName, ".swc")
        return load_swc(fileName)
    elseif endswith(fileName, ".swc.bin")
        return load_swc_bin(fileName)
    elseif endswith(fileName, ".bin") 
        return load_bin(fileName) 
    else 
        error("only support .swc or .swc.bin file: $fileName")
    end
end

function save(self::Neuron, fileName::AbstractString)
    if endswith(fileName, ".swc")
        save_swc(self, fileName)
    elseif endswith(fileName, ".swc.bin")
        save_swc_bin(self, fileName)
    elseif endswith(fileName, ".bin")
        save_bin(self, fileName) 
    else 
        error("unsupported format: $(fileName)")
    end 
end  

function load_swc( fileName::AbstractString )
    swcString = read( fileName , String)
    Neuron( swcString )
end

function save_swc(self::Neuron, fileName::AbstractString)
    swc = SWC(self)
    SWCs.save( swc, fileName )
end 

function load_swc_bin( fileName::AbstractString )
    swc = SWCs.load_swc_bin( fileName )
    Neuron( swc )
end 

function save_swc_bin( self::Neuron, fileName::AbstractString )
    SWCs.save_swc_bin( SWCs.SWC(self), fileName )
end 

function load_bin(fileName::AbstractString)
    open(fileName) do f 
        return Serialization.deserialize(f)
    end 
end  

function save_bin(self::Neuron, fileName::AbstractString)
    open(fileName, "w") do f 
        Serialization.serialize(f, self) 
    end  
    nothing 
end  
####################### properties ##############
"""
    get_path_length_normalized_features(self::Neuron)
get feature vector, returned as NamedTuple
Notethat some features were normalized using total path length 
"""
function get_path_length_normalized_features(self::Neuron)
    totalPathLength = get_total_path_length(self)
    return (
     branchingAngle = mean(map(i->get_branching_angle(self,i), 1:length(self))),
     tortuosity = mean(map(Segments.get_tortuosity, self)),
     symmetry = get_asymmetry(self),
     typicalRadius = get_typical_radius(self),
     numBranchingPoints = get_num_branching_points(self) / totalPathLength,
     medianBranchPathLength = median(get_segment_path_length_list(self)),
     surfaceArea = get_surface_area(self) / totalPathLength,
     volume = get_volume(self) / totalPathLength 
    )
end 


"""
    get_root_segment_id( self::Neuron )

the first segment should be the root segment, and the first node should be the root node 
"""
@inline function get_root_segment_id(self::Neuron) 1 end 
@inline function get_root_segment(self::Neuron) 
    get_segment_list(self)[get_root_segment_id(self)]
end 

"""
    get_root_node( self::Neuron )
the returned root node is a tuple of Float64, which represents x,y,z,r
"""
function get_root_node( self::Neuron ) 
    rootSegmentId = get_root_segment_id(self)
    rootSegment = get_segment_list(self)[ rootSegmentId ]
    rootSegment[1]
end 

"""
    get_num_segments(self::Neuron)
"""
@inline function get_num_segments(self::Neuron) length(self.segmentList) end

"""
    get_num_branching_points(self::Neuron)
"""
function get_num_branching_points(self::Neuron)
    numSegmentingPoint = 0
    for index in 1:get_num_segments(self)
        childrenSegmentIdList = get_children_segment_id_list(self, index)
        if length(childrenSegmentIdList) > 0
            numSegmentingPoint += 1
        end 
    end 
    numSegmentingPoint
end 

@inline function get_segment_list(self::Neuron{T}) where T
    self.segmentList::Vector{Segment{T}} 
end 

@inline function get_connectivity_matrix(self::Neuron) 
    self.connectivityMatrix 
end 

"""
    get_segment_order_list( self::Neuron )
following the Julia indexing style, the root segment order is 1.
"""
function get_segment_order_list(self::Neuron)
    segmentOrderList = Vector{Int}()
    sizehint!(segmentOrderList, get_num_segments(self))
    
    index2order = Dict{Int,Int}()
    # the first one is segment index, the second one is segment order 
    seedSegmentIdOrderList = Vector{NTuple{2,Int}}([(1,1)])
    while !isempty( seedSegmentIdOrderList )
        segmentId, segmentOrder = pop!(seedSegmentIdOrderList)
        index2order[segmentId] = segmentOrder 
        childrenSegmentIdList = get_children_segment_id_list(self, segmentId)
        for childSegmentId in childrenSegmentIdList 
            push!(seedSegmentIdOrderList, (childSegmentId, segmentOrder+1))
        end 
    end 
    
    for index in 1:get_num_segments( self )
        push!(segmentOrderList, index2order[index])
    end 
    segmentOrderList
end 

"""
    get_segment_length_list( self::Neuron )

get a vector of Integer, which represent the length of each segment 
"""
@inline function get_segment_length_list( self::Neuron ) map(length, get_segment_list(self)) end

"""
    get_node_num( self::Neuron )
get total number of nodes  
"""
@inline function get_node_num( self::Neuron )
    sum(map(length, get_segment_list(self)))
end 

"""
    get_node_list(self::Neuron)
get the node list. the first one is root node.
"""
function get_node_list( self::Neuron{T} ) where T
    nodeList = Vector{NTuple{4, T}}(undef, get_node_num(self))
    
    i = 0
    for segment in get_segment_list(self)
        for node in Segments.get_node_list(segment)
            i += 1
            nodeList[i] = node 
        end
    end 
    nodeList
end 

"""
    get_node_class_list(self::Neuron)
the class defines the type of node, such as axon, dendrite and soma.
"""
function get_node_class_list( self::Neuron )
    nodeClassList = Vector{UInt8}(undef, get_node_num(self))
    
    i = 0
    for segment in get_segment_list(self)
        nodeNumInSegment = length(segment)
        class = Segments.get_class(segment)
        nodeClassList[i+1:i+nodeNumInSegment] .= class
        i += nodeNumInSegment
    end 
    nodeClassList
end 

"""
    get_edge_list( self::Neuron )
get the edges with type of Vector{NTuple{2, Int}}
"""
function get_edge_list( self::Neuron )
    edgeList = Vector{NTuple{2,Int}}()
    segmentStartNodeIdList = Vector{Int}()
    segmentStopNodeIdList = Vector{Int}()
    # total number of nodes
    nodeNum = 0
    for (segmentId, segment) in enumerate(get_segment_list(self))
        push!(segmentStartNodeIdList, nodeNum+1)
        # build the edges inside segment
        for nodeId in nodeNum+1:nodeNum+length(segment)-1
            push!(edgeList, (nodeId, nodeId+1))
        end 
        # update node number
        nodeNum += length(segment)
        push!(segmentStopNodeIdList,  nodeNum)
    end 
    # add segment connections
    parentSegmentIdList, childSegmentIdList, _ = findnz( get_connectivity_matrix(self) )
    for (index, parentSegmentId) in enumerate( parentSegmentIdList )
        childSegmentId = childSegmentIdList[ index ]
        parentNodeId = segmentStopNodeIdList[ parentSegmentId ]
        childNodeId  = segmentStartNodeIdList[ childSegmentId ]
        push!( edgeList, (parentNodeId, childNodeId) )
    end 
    edgeList 
end 

function get_edge_num(self::Neuron)
    length( get_edge_list(self) )
end 

function get_path_to_soma_length(self::Neuron{T}, synapse::Synapse{T}) where T
    mergingSegmentId, closestNodeId = find_closest_node( self, synapse )
    get_path_to_soma_length( self, mergingSegmentId; nodeId=closestNodeId )
end 

"""
    get_path_to_soma_length(self::Neuron{T}, segmentId::Integer; 
                        segmentList::Vector{Segment}=get_segment_list(self), 
                        nodeId::Int = length(segmentList[segmentId]), 
                        segmentPathLengthList::Vector{T} = get_segment_path_length_list(self)) where T
"""
function get_path_to_soma_length(self::Neuron{T}, segmentId::Integer; 
                        segmentList::Vector{Segment{T}}=get_segment_list(self), 
                        nodeId::Integer = length(segmentList[segmentId]), 
                        segmentPathLengthList::Vector{T} = get_segment_path_length_list(self)) where T

    path2RootLength = Segments.get_path_length( segmentList[segmentId]; nodeId=nodeId )
    while true 
        parentSegmentId = get_parent_segment_id(self, segmentId )
        if parentSegmentId < 1 
            # root segment do not have parent 
            break 
        else
            path2RootLength += segmentPathLengthList[ parentSegmentId ]
            segmentId = parentSegmentId 
        end 
    end
    path2RootLength 
end

get_path_to_root_length = get_path_to_soma_length 

function get_pre_synapse_to_soma_path_length_list(self::Neuron{T}; 
                        segmentList::Vector{Segment{T}}=get_segment_list(self),
                        segmentPathLengthList::Vector{T}=get_segment_path_length_list(self)) where T
    preSynapseToSomaPathLengthList = Vector{T}()
    for (segmentId, segment) in enumerate( segmentList )
        preSynapseList = Segments.get_pre_synapse_list( segment )

        # get the node IDs with synapses attached
        nodeIdList = Int[]
        for (i,synapses) in enumerate(preSynapseList)
            if !ismissing(synapses)
                push!(nodeIdList, i)
            end
        end

        pathToSomaLengthList = map(nodeId->get_path_to_soma_length(self, segmentId; 
                                    segmentList=segmentList,
                                    segmentPathLengthList=segmentPathLengthList,
                                    nodeId=nodeId), nodeIdList)
        append!(preSynapseToSomaPathLengthList, pathToSomaLengthList)
    end
    preSynapseToSomaPathLengthList 
end 

function get_post_synapse_to_soma_path_length_list(self::Neuron{T}; 
                        segmentList::Vector{Segment{T}}=get_segment_list(self),
                        segmentPathLengthList::Vector{T}=get_segment_path_length_list(self)) where T
    postSynapseToSomaPathLengthList = Vector{T}()
    for (segmentId, segment) in enumerate( segmentList )
        postSynapseList = Segments.get_post_synapse_list( segment )

        # get the node IDs with synapses attached
        nodeIdList = Int[]
        for (i,synapses) in enumerate(postSynapseList)
            if !ismissing(synapses)
                push!(nodeIdList, i)
            end
        end

        pathToSomaLengthList = map(nodeId->get_path_to_soma_length(self, segmentId; 
                                    segmentList=segmentList,
                                    segmentPathLengthList=segmentPathLengthList,
                                    nodeId=nodeId), nodeIdList)
        append!(postSynapseToSomaPathLengthList, pathToSomaLengthList)
    end
    postSynapseToSomaPathLengthList 
end 


"""
get_segment_path_length_list(self::Neuron{T}; 
                                segmentList::Vector{Segment{T}}=get_segment_list(self),
                                class::{Nothing,UInt8}=nothing)
get euclidean path length of each segment 
"""
function get_segment_path_length_list( self::Neuron{T}; 
                                      segmentList::Vector{Segment{T}} = get_segment_list(self),
                                      class::Union{Nothing,UInt8}=nothing ) where T
    ret = Vector{T}()
    for (index, segment) in enumerate( segmentList )
        segmentPathLength = Segments.get_path_length( segment )
        # add the edge length to parent node
        parentSegmentId = get_parent_segment_id(self, index)
        if parentSegmentId > 0
            parentSegment = segmentList[ parentSegmentId ]
            parentNode = parentSegment[ end ]
            node = segment[1]
            segmentPathLength += norm( [map((x,y)->x-y, node[1:3], parentNode[1:3])...] )
        end

        if class == nothing || class==Segments.get_class(segment)
            # include all segment the class matches
            push!(ret, segmentPathLength)
        end
    end 
    ret 
end 

function get_node_distance_list(neuron::Neuron{T}) where T
    nodeList = Neurons.get_node_list(neuron)
    edgeList = Neurons.get_edge_list(neuron)
    nodeDistanceList = Vector{T}()
    sizehint!(nodeDistanceList, length(edgeList))
    for (src, dst) in edgeList
        d = norm( [map(-, nodeList[src][1:3], nodeList[dst][1:3])...] )
        push!(nodeDistanceList, d)
    end 
    nodeDistanceList
end 

"""
    get_total_path_length(self::Neuron; class::Union{Nothing,UInt8}=nothing)
the default class=nothing will include all of the segments
"""
@inline function get_total_path_length(self::Neuron; class::Union{Nothing,UInt8}=nothing)
    get_segment_path_length_list(self; class=class) |> sum
end 

"""
    get_segment_end_node_id_list( self::Neuron )

get a vector of integer, which represent the node index of segment end 
"""
@inline function get_segment_end_node_id_list( self::Neuron ) cumsum( get_segment_length_list(self) ) end 

"""
    get_num_nodes(self::Neuron)

get number of nodes 
"""
@inline function get_num_nodes( self::Neuron )
    segmentList = get_segment_list(self)
    sum(map(length, segmentList))
end 

"""
    get_furthest_terminal_node_pair_direction(self::Neurons)

Return: 
    vec::Vector{T}, the 3d vector representing the direction 
    maxNodeDistance::T, the maximum distance between terminal nodes 
"""
@inline function get_furthest_terminal_node_pair_direction(self::Neuron{T}) where T
    terminalNodeList = Neurons.get_terminal_node_list(self)
    get_furthest_terminal_node_pair_direction(terminalNodeList)
end 

function get_furthest_terminal_node_pair_direction(terminalNodeList::Vector{NTuple{4,T}}) where T
    # find the farthest node pair
    vec = zeros(T, 2)
    maxNodeDistance = zero(T)
    for i in 1:length(terminalNodeList)
        node1 = terminalNodeList[i][1:3]
        for j in i+1:length(terminalNodeList)
            node2 = terminalNodeList[j][1:3]
            v = [node1...] .- [node2...]
            d = norm(v)
            if d > maxNodeDistance
                vec = v
                maxNodeDistance = d
            end
        end
    end
    return vec, maxNodeDistance 
end 


"""
    get_mass_center( self::Neuron )
mass center was computed simply as center of all nodes 
"""
function get_mass_center( self::Neuron )
    nodeList = get_node_list(self)
    x = sum(map(n->n[1], nodeList)) / length(nodeList)
    y = sum(map(n->n[2], nodeList)) / length(nodeList)
    z = sum(map(n->n[3], nodeList)) / length(nodeList)
    (x,y,z)
end

"""
    get_typical_radius( self::Neuron )
Typical radius is the root-mean-square distance of dendritic arbor points to the center of mass (in nm)
"""
function get_typical_radius( self::Neuron )
    massCenter = get_mass_center( self )
    nodeList = get_node_list( self )
    typicalRadius = 0.0
    for node in nodeList 
        typicalRadius += norm([massCenter...] .- [node[1:3]...]) / length(nodeList)
    end 
    typicalRadius 
end 

"""
    get_asymmetry( self::Neuron )
asymmetry was measured by the euclidean distance between root node and mass center 
"""
function get_asymmetry( self::Neuron )
    massCenter = get_mass_center( self )
    root = get_root_node( self )
    norm( [massCenter...] .- [root[1:3]...] )
end 

"""
    get_children_segment_id_list(self::Neuron, parentSegmentId::Integer)
"""
@inline function get_children_segment_id_list(self::Neuron, parentSegmentId::Integer)
    get_children_segment_id_list(self.connectivityMatrix, parentSegmentId)
end

@inline function get_children_segment_id_list(connectivityMatrix::SparseMatrixCSC, 
                                      parentSegmentId::Integer)
    childrenSegmentIdList,_ = findnz(connectivityMatrix[parentSegmentId, :])
    childrenSegmentIdList 
end 

@inline function get_parent_segment_id( self::Neuron, childSegmentId::Integer )
    get_parent_segment_id(self.connectivityMatrix, childSegmentId) 
end

function get_parent_segment_id( connectivityMatrix::SparseMatrixCSC, childSegmentId::Integer )

    parentSegmentIdList,_ = findnz(connectivityMatrix[:, childSegmentId])
    if isempty( parentSegmentIdList ) 
        # no parent, this is a root segment
        return 0
    else
        if length(parentSegmentIdList) > 1
            @show parentSegmentIdList
            error("the neuron is supposed to be a tree, can not have loop with multiple parent segments.")
        end
        return parentSegmentIdList[1]
    end 
end

"""
    get_subtree_segment_id_list( self, segmentInde )
get the segment index list of subtree 
"""
function get_subtree_segment_id_list( self::Neuron, segmentId::Integer )
    @assert segmentId > 0 && segmentId <= get_num_segments(self)
    subtreeSegmentIdList = Vector{Integer}()
    seedSegmentIdList = [segmentId]
    while !isempty( seedSegmentIdList )
        seedSegmentId = pop!( seedSegmentIdList )
        push!(subtreeSegmentIdList, seedSegmentId)
        childrenSegmentIdList = get_children_segment_id_list( self, seedSegmentId )
        append!(seedSegmentIdList, childrenSegmentIdList)
    end 
    subtreeSegmentIdList 
end 

function get_terminal_segment_id_list( self::Neuron; startSegmentId::Integer = 1 )
    terminalSegmentIdList = Vector{Int}()
    for segmentId in startSegmentId:length(get_segment_list(self))
        childrenSegmentIdList = get_children_segment_id_list(self, segmentId)
        if isempty(childrenSegmentIdList) 
            push!(terminalSegmentIdList, segmentId)
        end 
    end 
    terminalSegmentIdList 
end

@inline function get_terminal_node_list( self::Neuron; startSegmentId::Integer = 1 )
    terminalSegmentIdList = get_terminal_segment_id_list( self )
    map( x -> Segments.get_node_list(x)[end], get_segment_list(self)[ terminalSegmentIdList ] )
end 

"""
    get_branching_angle( self::Neuron, segmentId::Integer; nodeDistance::AbstractFloat  )
if the node is too close the angle might not be accurate. For example, if they are voxel neighbors, the angle will alwasys be 90 or 45 degree. Thus, we have nodeDistance here to help make the nodes farther away.
Note that the acos returns angle in the format of radiens.
"""
function get_branching_angle( self::Neuron, segmentId::Integer; nodeDistance::Real = 5000.0 )
    segment = self[segmentId]
    parentSegmentId = get_parent_segment_id(self, segmentId)
    if parentSegmentId < 1
        # this is root segment
        return 0.0
    end
    parentSegment = self[ get_parent_segment_id(self, segmentId) ]
    if length(parentSegment) == 1 || length(segment) == 1
        return 0.0
    end 
    branchingNode = parentSegment[end]
    parentNode = parentSegment[end-1]
    for node in parentSegment 
        if Segments.get_nodes_distance(node, branchingNode) < nodeDistance
            parentNode = node 
            break 
        end 
    end
    if parentNode == branchingNode
        @warn("parent node is the same with branching node: $(branchingNode)")
        return 0.0
    end 

    childNode = segment[1]
    for index in length(segment):1
        node = segment[index]
        if Segments.get_nodes_distance(node, branchingNode) < nodeDistance 
            childNode = node 
            break 
        end 
    end 
    if childNode == branchingNode 
        @warn("child node is the same with branching node: $(branchingNode)")
        return 0.0 
    end 

    # compute the angle among the three nodes using the definition of dot product
    # get the two vectors
    v1 = map(-, parentNode[1:3], branchingNode[1:3] )
    v2 = map(-, branchingNode[1:3], childNode[1:3] )
    # normalize the vector
    nv1 = normalize([v1...]) 
    nv2 = normalize([v2...])
    dotProduct = dot(nv1, nv2)
    # tolerate some numerical varition. the dot product could go greater than 1.
    @assert dotProduct < 1.001 "impossible dotProduct: $(dotProduct)"
    return acos( min(1.0, dotProduct) )
end 

"""
    get_sholl_number(self::Neuron, shollRadius::AbstractFloat)
suppose there is a sphere centered on the root node. The sholl number is the number of contact points of the neuron and sphere. 
"""
function get_sholl_number(self::Neuron, shollRadius::AbstractFloat)
    root = get_root_node(self)
    nodeList = get_node_list( self )
    edgeList = get_edge_list( self )
    get_sholl_number( root, nodeList, edgeList, shollRadius )
end 

function get_sholl_number(root::NTuple{4, Float32}, nodeList::Vector{NTuple{4,Float32}}, 
                          edgeList::Vector{NTuple{2,Int}}, shollRadius::AbstractFloat)
    shollNum = 0
    for edge in edgeList 
        node1 = nodeList[ edge[1] ]
        node2 = nodeList[ edge[2] ]
        if  (Segments.get_nodes_distance(root, node1) > shollRadius) != (Segments.get_nodes_distance(root, node2) > shollRadius) 
            shollNum += 1
        end 
    end 
    shollNum 
end 

function get_sholl_number_list(self::Neuron, shollRadiusList::Vector)
    root = get_root_node(self)
    nodeList = get_node_list( self )
    edgeList = get_edge_list( self )
    shollNumList = zeros(Int, length(shollRadiusList))
    for (index, shollRadius) in enumerate(shollRadiusList)
        shollNumList[index] = get_sholl_number(root, nodeList, edgeList, shollRadius) 
    end 
    shollNumList 
end

"""
    get_sholl_number_list(self::Neuron, radiusStep::AbstractFloat; numStep::Integer=0)
if the number of steps is 0, we'll auto compute the step number according to the farthest terminal node.
"""
function get_sholl_number_list(self::Neuron, radiusStep::Real;
                               numStep::Integer=0)
    shollNumList = Vector{Float64}()
    if numStep==0
        # automatically computing step number 
        root = get_root_node(self)
        terminalNodeList = get_terminal_node_list(self)
        terminal2RootDistanceList = map(n->Segments.get_nodes_distance(root,n), terminalNodeList)
        maxDistance = maximum( terminal2RootDistanceList )
        numStep = fld(maxDistance, radiusStep)
    end

    radius = 0.0
    for n in 1:numStep
        radius += n*radiusStep 
        shollNum = get_sholl_number( self, radius )
        push!(shollNumList, shollNum)
    end 
    shollNumList 
end 

"""
    get_gyration_radius( self::Neuron )
The radius of gyration of an object is the square root of the sum of the squares of the radii from the center of mass to all the points on the object, divided by the square root of the number of points. 
"""
function get_gyration_radius( self::Neuron;  nodeList = get_node_list(self), 
                                                massCenter = get_mass_center(self) )
    distance2MassCenterList = map( n -> norm([massCenter...].-[n[1:3]...]), nodeList )
    sqrt( sum( distance2MassCenterList.^2 ) ) / sqrt(length(nodeList))
end 

"""
    get_surface_area(self::Neuron)
frustum-based 
http://www.analyzemath.com/Geometry_calculators/surface_volume_frustum.html
"""
function get_surface_area(self::Neuron)
    ret = zero(Float32) 
    segmentList = get_segment_list(self) 
    for (segmentId, segment) in enumerate(segmentList)
        ret += Segments.get_surface_area(segment) 
        parentSegmentId = get_parent_segment_id(self, segmentId)
        if parentSegmentId  > 0
            parentNode = segmentList[parentSegmentId][end]
            h = Segments.euclidean_distance(parentNode[1:3], segment[1][1:3]) 
            # average diameter 
            r1 = parentNode[4] 
            r2 = segment[1][4]
            ret += pi * (r1+r2) * sqrt( h*h + (r1-r2)*(r1-r2) )
        end 
    end 
    ret 
end

"""
    get_volume(self::Neuron)
frustum based volume:
http://jwilson.coe.uga.edu/emt725/Frustum/Frustum.cone.html
"""
function get_volume(self::Neuron)
    ret = zero(Float32)
    segmentList = get_segment_list(self)
    for (segmentId, segment) in enumerate(segmentList)
        ret += Segments.get_volume(segment)
        parentSegmentId = get_parent_segment_id(
                                        self, segmentId)
        if parentSegmentId > 0
            parentNode = segmentList[parentSegmentId][end]
            node = segment[1]
            r1 = node[4]
            r2 = parentNode[4]
            h = Segments.euclidean_distance(node[1:3], parentNode[1:3])
            ret += pi * h * (r1*r1 + r1*r2 + r2*r2) / Float32(3)
        end 
    end
    ret 
end 

"""
    get_fractal_dimension( self::Neuron )
compute fractal dimension using cumulative-mass method.
https://www.sciencedirect.com/science/article/pii/016502709400115W

Return:
    fractalDimension::Int 
    radiusList::Float64 the scanning disk radii 
    averageMassList::Float64 the average mass inside the disk   
"""
function get_fractal_dimension( self::Neuron )
    nodeList = get_node_list( self )
    massCenter = get_mass_center(self)
    gyrationRadius = get_gyration_radius( self; nodeList=nodeList, massCenter=massCenter )
    
    # find the nodes inside the gyration radius
    diskCenterList = Vector{NTuple{4,Float32}}()
    for node in nodeList
        if Segments.get_nodes_distance( massCenter, node ) < gyrationRadius 
            push!(diskCenterList, node)
        end 
    end 
    
    # radius list of a sequence of concentric disks
    radiusNum = 10
    radiusStep = 0.5*gyrationRadius/radiusNum
    radiusList = zeros(radiusNum)
    radius = 0.0
    for i in 1:radiusNum
        radius += i*radiusStep 
        radiusList[i] = radius 
    end 

    # iterate all the nodes inside gyration radius as the center of scanning disks
    averageMassList = zeros( length(radiusList) )
    for (centerId, center) in enumerate(diskCenterList)
        for (radiusId, radius) in enumerate(radiusList)
            for node in nodeList 
                if Segments.get_nodes_distance( center,  node) < radius
                    averageMassList[radiusId] += 1.0 / length(diskCenterList)
                end 
            end 
        end 
    end 

    # fit the curve and get slop as fractal dimension
    model(x,p) = p[1]*x + p[2]
    p0 = [1.0, 0]
    fit = curve_fit(model, log(radiusList), log.(averageMassList), p0)
    fractalDimension = fit.param[1]
    return fractalDimension, radiusList, averageMassList 
end 


"""
	get_mask(self::Neuron, voxelSize::Union{Tuple, Vector}) 
compute binary representation of the neuron skeleton
"""
function get_mask(self::Neuron, voxelSize::Union{Tuple, Vector})
    nodeList = get_node_list(self)
    voxelSet = Set{NTuple{3, Int}}()
    for node in nodeList 
        voxelCoordinate = map((x,y)->round(Int, x/y), node[1:3], voxelSize)
        push!(voxelSet, (voxelCoordinate...,))
    end 
    boundingBox = Segments.BoundingBox( voxelSet )
    range = BoundingBoxes.get_unit_range(boundingBox)
    sz = size(boundingBox)
    # initialize the map 
    mask = zeros(Bool, sz)
    mask = OffsetArray(mask, range...)

    @inbounds for voxel in voxelSet
        # add a small number to make sure that there is no 0
        mask[voxel...] = true 
    end 
    mask
end 

function get_bounding_box(self::Neuron)
    segmentList = get_segment_list(self)
    bb = Segments.get_bounding_box( segmentList[1] )
    for segment in segmentList[2:end]
        bb = union(bb, Segments.get_bounding_box(segment))
    end 
    return bb 
end 

function get_2d_binary_projection(self::Neuron; axis=3, voxelSize=VOXEL_SIZE)
    nodeList = get_node_list(self)
    voxelSet = Set{NTuple{3, Int}}()
    for node in nodeList 
        voxelCoordinate = map((x,y)->round(Int, x/y), node[1:3], voxelSize)
        push!(voxelSet, (voxelCoordinate...,))
    end 
    boundingBox = Segments.BoundingBox( voxelSet ) 
    range = Segments.BoundingBoxes.get_unit_range(boundingBox)[1:2]
    sz = size(boundingBox)[1:2]
    # initialize the map                                                   
    mask = zeros(Bool, sz)
    mask = OffsetArray(mask, range...)

    @inbounds for voxel in voxelSet 
        # add a small number to make sure that there is no 0 
        mask[voxel[1:2]...] = true 
    end 
    mask
end 

"""
    get_arbor_density_map(self::Neuron, 
                        voxelSize::Union{Tuple, Vector},
                        gaussianFilterStd ::AbstractFloat)
compute the arbor density map

Return:
    * densityMap::OffsetArray 
"""
function get_arbor_density_map(self::Neuron, voxelSize::Union{Tuple, Vector},
                                            gaussianFilterStd::AbstractFloat)
    densityMap = get_mask(self, voxelSize)
    # gaussion convolution
    kernel = Kernel.gaussian((3,3,3))
    # kernel = fill(gaussianFilterStd, 3) |> Kernel.gaussian
    println("convolution of gaussian kernel...")
    @time densityMap = imfilter(densityMap, kernel)
    # normalize by total path length of the neuron
    #totalPathLength = get_total_path_length(self)
    #densityMap .*= totalPathLength / norm(densityMap[:]) 
    # densityMap ./= norm(densityMap[:]) 
    densityMap
end

"""
    translate_soma_to_coordinate_origin(self::Neuron, densityMap::OffsetArray)
translate the soma to coordinate (0,0,0)
"""
function translate_soma_to_coordinate_origin(self::Neuron, densityMap::OffsetArray, 
                                             voxelSize::Tuple)
    # assume that the first node is the soma node with largest radius
    somaCorrdinate = get_root_node(self)[1:3]
    newRange = map((x,y,s)->(x .- round(Int,y/s)), 
                   axes(densityMap), somaCorrdinate, voxelSize)
    OffsetArray(parent(densityMap), newRange...)
end 

"""
    get_arbor_density_map_distance(self::OffsetArray, other::OffsetArray)
compute the arbor density map distance. Note that the two density maps should both be 
translated to make the soma located in (0,0,0).
"""
function get_arbor_density_map_distance(self::OffsetArray{T,N,Array{T,N}}, 
                                        other::OffsetArray{T,N,Array{T,N}}) where {T,N}
    # find the intersection of two offset arrays
    intersectIndices = map(intersect, axes(self), axes(other))
    d = 0.0
    # add the distance of intersection
    d += sum(abs2.(self[intersectIndices...] .- other[intersectIndices...])) 
    # add the distance of self from non-intersect regions
    selfDiff  = map(setdiff, axes(self),  intersectIndices)
    otherDiff = map(setdiff, axes(other), intersectIndices)
    d += sum(abs2.(self[selfDiff...]))
    d += sum(abs2.(other[otherDiff...]))
    d
end 

function _find_correlation_peak_along_axis(self::Array{T,3}, other::Array{T,3}, 
                                                                    axis::Int) where {T}
	crossCorrelation = _cross_correlation_along_axis(self, other, axis)
   	maximum(crossCorrelation)    
end 

"""                                                                                                                                            
    get_arbor_density_map_overlap_min_distance(self ::OffsetArray, other::OffsetArray)                                                         
                                                                                                                                               
translate the density map and find the place with minimum distance in overlapping region.                                                      
use cross correlation in 2D to find the position with maximum value as approximated location.                                                  
""" 
function get_arbor_density_map_overlap_min_distance(self ::Array{T,3},
                                                    other::Array{T,3}) where {T}
    crossCorrelationZ = _find_correlation_peak_along_axis(self, other, 3)
    crossCorrelationY = _find_correlation_peak_along_axis(self, other, 2)
    crossCorrelationX = _find_correlation_peak_along_axis(self, other, 1)
    one(typeof(crossCorrelationX)) - mean([crossCorrelationX, crossCorrelationY, crossCorrelationZ]) 
end 

############################### Base functions ###################################

function Base.show(self::Neuron)
    println("neuron with ", length(self), " segments, ", get_node_num(self), " nodes, ", 
            get_num_pre_synapses(self), " buttons, and ",
            get_num_post_synapses(self), " post synapses.")
    nothing
end 

function Base.iterate(self::Neuron, state::Int=1)
    if state > length(self)
        return nothing 
    end 

    return self[state], state+1
end 

@inline function Base.length(self::Neuron)
    length(get_segment_list(self))
end 

@inline function Base.getindex(self::Neuron, index::Integer)
    get_segment_list(self)[index]
end

"""
    Base.merge(self::Neuron, other::Neuron, mergingSegmentId::Integer, mergingNodeIdInSegment::Integer)

merge two nets at a specific segment location.
the root node of second net will be connected to the first net
"""
function Base.merge(self::Neuron{T}, other::Neuron{T}, 
                    mergingSegmentId::Integer, mergingNodeIdInSegment::Integer) where T
    @assert !isempty(self)
    @assert !isempty(other)
    @assert mergingSegmentId > 0
    segmentList1 = get_segment_list( self  )
    segmentList2 = get_segment_list( other )
    num_segments1 = get_num_segments(self)
    num_segments2 = get_num_segments(other)
    
    if mergingNodeIdInSegment == length(segmentList1[mergingSegmentId])
        println("connecting to a segment end...")
        childrenSegmentIdList = get_children_segment_id_list(self, mergingSegmentId)
        if length(childrenSegmentIdList) > 0
            println("the closest segment have children, do not merge the root segment")
            total_num_segments = num_segments1 + num_segments2 
            mergedSegmentList = vcat(segmentList1, segmentList2)
            @assert length(mergedSegmentList) == total_num_segments 
            mergedConnectivityMatrix = 
                        spzeros(Bool, total_num_segments, total_num_segments)
            mergedConnectivityMatrix[   1:size(self.connectivityMatrix,1), 
                                        1:size(self.connectivityMatrix,2)] = 
                                                                self.connectivityMatrix 
            # do not include the connection of root in net2
            mergedConnectivityMatrix[num_segments1+1 : end, 
                                     num_segments1+1 : end] = other.connectivityMatrix 
            mergedConnectivityMatrix[mergingSegmentId, num_segments1+1] = true
            return Neuron{T}(mergedSegmentList, mergedConnectivityMatrix)
        else 
            println("the closest segment is a terminal segment, merge the root segment of the other net")
            if num_segments2 == 1
                #  only one segment of the other net
                mergedSegmentList = copy(segmentList1)
                mergedSegmentList[mergingSegmentId] = 
                            merge(mergedSegmentList[mergingSegmentId], segmentList2[1])
                            return Neuron{T}(mergedSegmentList, self.connectivityMatrix)
            else 
                # the connection point is the end of a segment, no need to break segment
                # merge the root segment of second net to the first net
                # assume that the first segment is the root segment
                mergedSegmentList = vcat( segmentList1, segmentList2[2:end] )
                mergedSegmentList[mergingSegmentId] = 
                            merge(mergedSegmentList[mergingSegmentId], segmentList2[1])

                # total number of segments
                total_num_segments = num_segments1 + num_segments2 - 1 
                @assert length(mergedSegmentList) == total_num_segments 
                
                mergedConnectivityMatrix = 
                                spzeros(Bool, total_num_segments, total_num_segments)
                mergedConnectivityMatrix[
                                    1:size(self.connectivityMatrix,1), 
                                    1:size(self.connectivityMatrix,2)] = 
                                                        self.connectivityMatrix 
                # do not include the connection of root in net2
                mergedConnectivityMatrix[num_segments1+1 : end, 
                                         num_segments1+1 : end] = 
                                                other.connectivityMatrix[2:end, 2:end]
                # reestablish the connection of root2
                childrenSegmentIdList2 = get_children_segment_id_list(other, 1)
                for childSegmentId2 in childrenSegmentIdList2
                    mergedConnectivityMatrix[ mergingSegmentId, 
                                              num_segments1+childSegmentId2-1 ] = true 
                end
                return Neuron{T}(mergedSegmentList, mergedConnectivityMatrix)
            end 
        end 
    else 
        println("need to break the segment and rebuild connectivity matrix")
        total_num_segments = num_segments1 + 1 + num_segments2 
        mergedSegmentList = segmentList1 
        mergedConnectivityMatrix = spzeros(Bool, total_num_segments, total_num_segments)
        
        # need to break the segment and then stitch the new segments
        segmentPart1, segmentPart2 = split(segmentList1[mergingSegmentId], mergingNodeIdInSegment)

        mergedSegmentList[mergingSegmentId] = segmentPart1
        # the connectivity is the same with old net
        mergedConnectivityMatrix[1:size(self.connectivityMatrix,1), 
                                 1:size(self.connectivityMatrix,2)] = self.connectivityMatrix 

        # reconnect the breaked two segments
        push!(mergedSegmentList, segmentPart2)
        mergedConnectivityMatrix[mergingSegmentId, num_segments1+1] = true 

        # redirect the children segments to segmentPart2
        childrenSegmentIdList = get_children_segment_id_list(self, mergingSegmentId)
        for childSegmentId in childrenSegmentIdList
            # remove old connection
            mergedConnectivityMatrix[mergingSegmentId, childSegmentId] = false
            # build new connection
            mergedConnectivityMatrix[num_segments1+1, childSegmentId] = true 
        end 

        # merge the other net
        mergedSegmentList = vcat(mergedSegmentList, segmentList2)
        mergedConnectivityMatrix[
            num_segments1+2 : num_segments1+2+size(other.connectivityMatrix,1)-1, 
            num_segments1+2 : num_segments1+2+size(other.connectivityMatrix,2)-1] = other.connectivityMatrix

        # establish the connection between two nets
        mergedConnectivityMatrix[mergingSegmentId, num_segments1+2] = true

        # remove the false elements 
        dropzeros!(mergedConnectivityMatrix) 

        # assert the parent number to ensure that this is a tree structure without loop
        for i in 1:size(mergedConnectivityMatrix,1)
            parentSegmentIdList, _ = findnz(mergedConnectivityMatrix[:,i])
            if length(parentSegmentIdList) > 1 
                @show parentSegmentIdList 
            end 
            @assert length(parentSegmentIdList) <= 1
        end 

        # create new merged net
        return Neuron{T}(mergedSegmentList, mergedConnectivityMatrix)
    end
end 

function Base.isempty(self::Neuron)    isempty(self.segmentList) end 

"""
    Base.split(self::Neuron, splitSegmentId::Integer; nodeIdInSegment::Integer=0)

split the segment net into two neurons
the nodeIdInSegment will be included in the first main net including the original root node  
"""
function Base.split(self::Neuron{T}, splitSegmentId::Integer; nodeIdInSegment::Integer=1) where T
    if splitSegmentId == 1 && nodeIdInSegment < 1 
        return Neuron(), self
    end 

    segmentList = get_segment_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    segmentIdListInSubtree2 = get_subtree_segment_id_list( self, splitSegmentId )
    subtree1Segment, subtree2RootSegment = split(segmentList[splitSegmentId], 
                                                                nodeIdInSegment)

    # build the split out subtree2, which do not contain the root node
    segmentList2 = [subtree2RootSegment]
    append!(segmentList2, segmentList[ segmentIdListInSubtree2 ])
    segmentIdMap = Dict{Int,Int}()
    sizehint!(segmentIdMap, length(segmentIdListInSubtree2))
    for (segmentIdInSubtree2, mainTreeSegmentId ) in 
                                                    enumerate(segmentIdListInSubtree2)
        segmentIdMap[ mainTreeSegmentId ] = segmentIdInSubtree2 
    end
    connectivityMatrix2 = spzeros( Bool, length(segmentList2), length(segmentList2) )
    for (segmentIdInSubtree2, segmentIdInOriginalTree) in 
                                            enumerate(segmentIdListInSubtree2)
        parentSegmentIdInOriginalTree = 
                                get_parent_segment_id(self, segmentIdInOriginalTree)
        if haskey(segmentIdMap, parentSegmentIdInOriginalTree)
            parentSegmentIdInSubtree2 = segmentIdMap[ parentSegmentIdInOriginalTree ]
            connectivityMatrix2[ parentSegmentIdInSubtree2, segmentIdInSubtree2 ] = true 
        else 
            # this is the root segment of new tree 
            println("find root of new tree in the original tree : $(parentSegmentIdInOriginalTree)")
        end 
    end 
    subtree2 = Neuron{T}( segmentList2, connectivityMatrix2 )

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    segmentList1 = Vector{Segment{T}}()
    segmentIdListInSubtree1 = Vector{Int}()
    segmentIdMap = Dict{Int,Int}()
    # the first one is the root segment
    seedSegmentIdList = [1]
    # build segmentList1 without the splitout segment, will deal with it later
    while !isempty( seedSegmentIdList )
        seedSegmentId = pop!( seedSegmentIdList )
        if seedSegmentId != splitSegmentId 
            push!(segmentList1, segmentList[seedSegmentId])
            push!(segmentIdListInSubtree1, seedSegmentId)
            # map the old segment index to the new segment index
            # 3 => 2 means the 3rd segment of old tree is the 2nd segment of new tree
            segmentIdMap[ seedSegmentId ] = length(segmentList1)
            childrenSegmentIdList = get_children_segment_id_list(self, seedSegmentId)
            append!(seedSegmentIdList, childrenSegmentIdList)
        else 
            # deal with the split out segment
            if nodeIdInSegment != 1
                @assert length(subtree1Segment) != 0
                push!(segmentList1, subtree1Segment)
                segmentIdMap[ seedSegmentId ] = length(segmentList1)
                # do not need to travase the subtree2 
                # so don't put children to seed list
            end 
        end 
    end 
    connectivityMatrix1 = spzeros(Bool, length(segmentList1), length(segmentList1))
    for (segmentIdInSubtree1, segmentIdInOriginalTree) in 
                                                    enumerate(segmentIdListInSubtree1)
        parentSegmentIdInOriginalTree = 
                                get_parent_segment_id(self, segmentIdInOriginalTree)
        if parentSegmentIdInOriginalTree > 0
            parentSegmentIdInSubtree1 = segmentIdMap[ parentSegmentIdInOriginalTree ]
            connectivityMatrix1[ parentSegmentIdInSubtree1, segmentIdInSubtree1 ] = true 
        end 
    end 
    subtree1 = Neuron{T}( segmentList1, connectivityMatrix1 )
    return subtree1, subtree2
end


function remove_segments(self::Neuron, removeSegmentIdList::Vector{Int})
    remove_segments(self, Set{Int}( removeSegmentIdList ))
end 

function remove_segments(self::Neuron{T}, removeSegmentIdSet::Set{Int}) where T
    if isempty(removeSegmentIdSet)
        return self 
    end 
    @assert !(1 in removeSegmentIdSet) "should not contain the root segment!"

    segmentList = get_segment_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    newSegmentList = Vector{Segment{T}}()
    newSegmentIdList = Vector{Int}()
    segmentIdMap = Dict{Int,Int}()
    # the first one is the root segment
    seedSegmentIdList = [get_root_segment_id(self)]
    # build segmentList1 without the splitout segment, will deal with it later
    while !isempty( seedSegmentIdList )
        seedSegmentId = pop!( seedSegmentIdList )
        if !(seedSegmentId in removeSegmentIdSet) 
            # should contain this segment
            push!(newSegmentList, segmentList[seedSegmentId])
            push!(newSegmentIdList, seedSegmentId)
            # map the old segment index to the new segment index
            segmentIdMap[ seedSegmentId ] = length(newSegmentList)
            childrenSegmentIdList = get_children_segment_id_list(self, seedSegmentId)
            append!(seedSegmentIdList, childrenSegmentIdList)
        end 
    end 
    newConnectivityMatrix = spzeros(Bool, length(newSegmentList), length(newSegmentList))
    for (newSegmentId, segmentIdInOriginalTree) in enumerate(newSegmentIdList)
        parentSegmentIdInOriginalTree = 
                                get_parent_segment_id(self, segmentIdInOriginalTree)
        if parentSegmentIdInOriginalTree > 0
            newParentSegmentId = segmentIdMap[ parentSegmentIdInOriginalTree ]
            newConnectivityMatrix[ newParentSegmentId, newSegmentId ] = true 
        end 
    end
    neuron = Neuron{T}( newSegmentList, newConnectivityMatrix )
    neuron = merge_only_child(neuron)
    return neuron
end 

"""
    merge_only_child(self::Neuron)
If a segment only have one child, there should not be a branching point here. 
It should be merged with it's only-child.
Find the segment mapping by depth-first search. 
"""
function merge_only_child(self::Neuron{T}) where T
    segmentList = get_segment_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    newSegmentList = Vector{Segment{T}}()
    # root always map to root 
    oldId2NewId = Dict{Int,Int}(0=>0)
    seedSegmentIdList = [1]
   
    # depth first traverse 
    while !isempty(seedSegmentIdList)
        seedSegmentId = pop!(seedSegmentIdList)
        seedSegment = segmentList[seedSegmentId]
        childrenSegmentIdList = get_children_segment_id_list(self, seedSegmentId)
        oldSegmentIdList = [seedSegmentId]    
        # find all the single child segments along this line 
        while length(childrenSegmentIdList) == 1
            childSegmentId = childrenSegmentIdList[1]
            push!(oldSegmentIdList, childSegmentId)
            childSegment = segmentList[childSegmentId]
            seedSegment = merge(seedSegment, childSegment)
            @assert length(seedSegment) > length(childSegment)
            childrenSegmentIdList = get_children_segment_id_list(self, 
                                                                 childSegmentId)
        end 
        @assert length(seedSegment) > 0
        push!(newSegmentList, seedSegment)
        for oldSegmentId in oldSegmentIdList 
            oldId2NewId[oldSegmentId] = length(newSegmentList)
        end 

        @assert length(childrenSegmentIdList) != 1
        for childSegmentId in childrenSegmentIdList
            push!(seedSegmentIdList, childSegmentId) 
        end
    end 
    
    # create new connectivity matrix 
    newConnectivityMatrix = spzeros(Bool, length(newSegmentList), length(newSegmentList))
    for oldSegmentId in 1:length(segmentList)
        oldParentSegmentId = get_parent_segment_id(connectivityMatrix, oldSegmentId)
        newParentSegmentId = oldId2NewId[oldParentSegmentId] 
        newSegmentId = oldId2NewId[oldSegmentId]
        if newParentSegmentId!=newSegmentId && newParentSegmentId>0 
            newConnectivityMatrix[newParentSegmentId, newSegmentId] = true 
        end 
    end
    Neuron{T}(newSegmentList, newConnectivityMatrix)
end 

"""
    remove_subtree_in_soma(self::Neuron)

remove the subtree which is inside the soma, which is an artifact of TEASAR algorithm
"""
function remove_subtree_in_soma( self::Neuron )
    removeSegmentIdList = Vector{Int}()
    rootSegmentId = get_root_segment_id( self )
    rootNode = get_root_node( self )
    childrenSegmentIdList = get_children_segment_id_list( self, rootSegmentId )
    for segmentId in childrenSegmentIdList 
        terminalNodeList = get_terminal_node_list( self; startSegmentId=segmentId )
        if all( map(n->Segments.get_nodes_distance(n, rootNode) < rootNode[4]*2, 
                                                                terminalNodeList ) )
            # println("remove segment: $segmentId")
            push!(removeSegmentIdList, segmentId)
        end 
    end
    return remove_segments( self, removeSegmentIdList )
end

function remove_hair( self::Neuron; radiiScale::Float32 = Float32(2),
                                    minDistance::Float32 = Float32(2000) )
    segmentList = get_segment_list(self)
    removeSegmentIdList = Vector{Int}()
    for terminalSegmentId in get_terminal_segment_id_list(self)
        parentSegmentId = get_parent_segment_id(self, terminalSegmentId)
        if parentSegmentId > 0
            # this is not a root segment
            parentSegment =segmentList[parentSegmentId] 
            parentNode = parentSegment[end]
            terminalNode = segmentList[ terminalSegmentId ][end]
            distance = Segments.get_nodes_distance(terminalNode, parentNode)
            if distance < radiiScale*parentNode[4] || distance < minDistance
                # println("remove segment $(terminalSegmentId) with distance of $(distance) and radius of $(parentNode[4])")
                push!(removeSegmentIdList, terminalSegmentId)
            end
        end 
    end 
    return remove_segments(self, removeSegmentIdList)
end

"""
    remove_terminal_blobs( self::Neuron )
some terminal segmentation was fragmented. a lot of blobs was attatched to dendrite.
The blob terminal segment path length is normally smaller than the distance to parent dendrite.
"""
function remove_terminal_blobs( self::Neuron )
    terminalSegmentIdList = get_terminal_segment_id_list( self )
    blobTerminalSegmentIdList = Vector{Int}()
    for index in terminalSegmentIdList 
        segment = self[index]
        distance2Parent = 0.0
        try 
            parentSegment = self[ get_parent_segment_id(self, index) ]
            distance2Parent = Segments.get_nodes_distance( segment[1], parentSegment[end] )
        catch err 
            @warn("this terminal segment is root segment!")
            @assert get_parent_segment_id(self, index) < 1
        end 
        segmentInnerPathLength = Segments.get_path_length( segment )
        if distance2Parent > segmentInnerPathLength
            # println("remove terminal blob. distance to parent: $(distance2Parent). segment inner path length: $(segmentInnerPathLength).")
            push!(blobTerminalSegmentIdList, index)
        end 
    end 
    return remove_segments( self, blobTerminalSegmentIdList )
end

"""
    remove_redundent_nodes( self::Neuron )
if neighboring node is the same, remove one of them 
"""
function remove_redundent_nodes(self::Neuron{T}) where T
    newSegmentList = Vector{Segment{T}}()
    sizehint!(newSegmentList, get_num_segments(self))
    removeSegmentIdList = Vector{Int}()
    segmentList = get_segment_list(self)
    for (segmentId, segment) in enumerate( segmentList )
        Segments.remove_redundent_nodes!( segment )
        parentSegmentId = get_parent_segment_id(self, segmentId)
        if parentSegmentId > 0 
            parentSegment = segmentList[ parentSegmentId ]
            if segment[1] == parentSegment[end]
                segment = Segments.remove_node(segment, 1)
            end
        end 
        push!(newSegmentList, segment) 
        if length(segment) == 0
            # empty segment, should be removed later
            push!(removeSegmentIdList, segmentId)
        end 
    end 
    newNeuron = Neuron{T}( newSegmentList, get_connectivity_matrix(self) )
    return remove_segments(newNeuron, removeSegmentIdList)
end


########################## type convertion ####################
"""
    NodeNets.NodeNet( self::Neuron )
transform to NodeNet, the first node is the root node.
"""
function NodeNets.NodeNet(self::Neuron)
    nodeList = get_node_list( self )
    edges = get_edge_list( self )
    nodeClassList = get_node_class_list(self)
    
    I = map(x->UInt32(x[1]), edges)
    J = map(x->UInt32(x[2]), edges)
    
    # the connectivity matrix should be symmetric
    connectivityMatrix = sparse([I..., J...,], [J..., I...,],true, 
                                length(nodeList), length(nodeList))
    NodeNet(nodeList, nodeClassList, connectivityMatrix)    
end 

function SWCs.SWC(self::Neuron)
    # initialize swc, which is a list of point objects
    swc = SWCs.SWC()
    # the node index of each segment, will be used for connecting the child segment
    segmentEndNodeIdList = get_segment_end_node_id_list(self)
    # connectivity matrix of segments
    segmentConnectivityMatrix = get_connectivity_matrix(self)

    for (segmentId, segment) in enumerate( get_segment_list(self) )
        for (nodeId, node) in enumerate( Segments.get_node_list( segment ))
            # type, x, y, z, r, parent
            # in default, this was assumed that the connection is inside a segment, 
            # the parent is simply the previous node 
            pointObj = SWCs.PointObj( Segments.get_class(segment), 
                                     node[1], node[2], node[3], node[4], length(swc) )
            
            if nodeId == 1
                # the first node should connect to other segment or be root node
                if length(swc) == 0
                    pointObj.parent = -1
                else
                    # find the connected parent segment index
                    parentSegmentId = get_parent_segment_id(self, segmentId)
                    # the node index of parent segment end
                    pointObj.parent= segmentEndNodeIdList[ parentSegmentId ]
                end 
            end 
            push!(swc, pointObj)
        end 
    end
    swc
end 

"""
get binary buffer formatted as neuroglancer nodeNet.

# Binary format
    UInt32: number of vertex
    UInt32: number of edges
    Array{Float32,2}: Nx3 array, xyz coordinates of vertex
    Array{UInt32,2}: Mx2 arrray, node index pair of edges
reference: 
https://github.com/seung-lab/neuroglancer/wiki/Skeletons
"""
function get_neuroglancer_precomputed(self::Neuron)
    node_list = get_node_list(self)
    edge_list = get_edge_list(self)

    # total number of bytes
    num_bytes = 4 + 4 + 4*3*length(node_list) + 4*2*length(edge_list)
    buffer = IOBuffer( read=false, write=true, maxsize=num_bytes )

    # write the number of vertex
    write(buffer, UInt32(length(node_list)))
    write(buffer, UInt32(length(edge_list)))
    # write the node coordinates
    for node in node_list
        write(buffer, node[1])
        write(buffer, node[2])
        write(buffer, node[3])
    end 
    for edge in edge_list
        write(buffer, UInt32( edge[1] - 1 ))
        write(buffer, UInt32( edge[2] - 1 ))
    end
    bin = Vector{UInt8}(take!(buffer))
    close(buffer)
    return bin 
end

"""
return all the pre synapses in a list 
"""
function get_all_pre_synapse_list(self::Neuron{T}) where T
    ret = Synapse{T}[]
    for segment in self 
        preSynapseList = Segments.get_pre_synapse_list(segment)
        for synapse in preSynapseList
            if !ismissing(synapse)
                append!(ret, synapse)
            end
        end
    end
    ret 
end 

"""
return all the post synapses in a list 
"""
function get_all_post_synapse_list(self::Neuron{T}) where T
    ret = Synapse{T}[]
    for segment in self 
        postSynapseList = Segments.get_post_synapse_list(segment)
        for synapse in postSynapseList
            if !ismissing(synapse)
                append!(ret, synapse)
            end
        end
    end
    ret 
end 

function get_num_pre_synapses(self::Neuron)
    ret = 0
    for segment in get_segment_list(self)
        ret += Segments.get_num_pre_synapses(segment)
    end 
    ret 
end 

function get_num_post_synapses(self::Neuron)
    ret = 0
    for segment in get_segment_list(self)
        ret += Segments.get_num_post_synapses(segment)
    end 
    ret 
end 

############################### manipulations ##########################################

"""
    attach_pre_synapses!(self::Neuron{T}, synapseList::Vector{Synapse{T}})

"""
function attach_pre_synapses!(self::Neuron{T}, synapseList::Vector{Synapse{T}}) where T
    kdTree, nodePosition = build_tree_with_position(self)
    for synapse in synapseList 
        attach_pre_synapse!(self, synapse, kdTree, nodePosition)
    end 

    for segment in self 
        Segments.adjust_class!(segment)
    end

    numSynapse = get_num_pre_synapses(self)
    if numSynapse < length(synapseList)
        @warn("we got $(length(synapseList)) synapses, but only $numSynapse were attached.")
    end 
    nothing 
end

"""
    attach_post_synapses!(self::Neuron{T}, synapseList::Vector{Synapse{T}})
"""
function attach_post_synapses!(self::Neuron{T}, synapseList::Vector{Synapse{T}}) where T
    kdTree, nodePosition = build_tree_with_position(self)
    for synapse in synapseList 
        attach_post_synapse!(self, synapse, kdTree, nodePosition)
    end 

    for segment in self 
        Segments.adjust_class!(segment)
    end 

    numSynapse = get_num_post_synapses(self)
    if numSynapse < length(synapseList)
        @warn("we got $(length(synapseList)) synapses, but only $numSynapse were attached.")
    end 
    nothing 
end 

"""
    attach_pre_synapses!(self::Neuron, synapseTable::DataFrame)
make sure that the synapse table contains *only* the presynapses of this neuron 
"""
function attach_pre_synapses!(self::Neuron, synapseTable::DataFrame)
    kdTree, nodePosition = build_tree_with_position(self)
    for row in DataFrames.eachrow( synapseTable )
        synapse = Synapse(row)
        attach_pre_synapse!(self, synapse, kdTree, nodePosition)
    end 
    # adjust segment class according to this new information
    for segment in get_segment_list(self)
        Segments.adjust_class!(segment)
    end 
    
    numSynapse1 = DataFrames.nrow(synapseTable)
    numSynapse2 = get_num_pre_synapses(self)
    if numSynapse2 < numSynapse1
        @warn("we got $(numSynapse1) synapses, but only $numSynapse2 were attached.")
    end 

    nothing 
end 

"""
    attach_post_synapses!(self::Neuron, synapseTable::DataFrame)
make sure that the synapse table contains *only* the postsynapses of this neuron 
"""
function attach_post_synapses!(self::Neuron, synapseTable::DataFrame)
    kdTree, nodePosition = build_tree_with_position(self)
    for row in DataFrames.eachrow( synapseTable )
        synapse = Synapse(row)
        attach_post_synapse!(self, synapse, kdTree, nodePosition)
    end  
    # adjust segment class according to this new information
    for segment in get_segment_list(self)
        Segments.adjust_class!(segment)
    end 
    
    numSynapse1 = DataFrames.nrow(synapseTable)
    numSynapse2 = get_num_post_synapses(self)
    if numSynapse2 < numSynapse1 
        @warn("we got $numSynapse1 synapses, but only $numSynapse2 were attached.")
    end 

    nothing
end 

function attach_pre_synapse!(self::Neuron, synapse::Segments.Synapse, 
                             kdTree::KDTree, nodePosition::Matrix{Int})
    synapticCoordinate = Synapses.get_pre_synaptic_coordinate( synapse )
    segmentId, nodeIdInSegment,_ = find_closest_node( kdTree, nodePosition, 
                                                        synapticCoordinate )
    segmentList = get_segment_list(self)
    Segments.attach_pre_synapse!(segmentList[segmentId], nodeIdInSegment, synapse)
    nothing
end 

function attach_pre_synapse!(self::Neuron, synapse::Segments.Synapse)
    synapticCoordinate = Synapses.get_pre_synaptic_coordinate( synapse )
    segmentId, nodeIdInSegment, _ = find_closest_node( self, synapticCoordinate )
    segmentList = get_segment_list(self)
    Segments.attach_pre_synapse!(segmentList[segmentId], nodeIdInSegment, synapse)
    nothing 
end 

function attach_post_synapse!(self::Neuron, synapse::Segments.Synapse, 
                             kdTree::KDTree, nodePosition::Matrix{Int})
    synapticCoordinate = Synapses.get_post_synaptic_coordinate( synapse )
    segmentId, nodeIdInSegment,_ = find_closest_node( kdTree, nodePosition, 
                                                        synapticCoordinate )
    segmentList = get_segment_list(self)
    Segments.attach_post_synapse!(segmentList[segmentId], nodeIdInSegment, synapse)
    nothing
end 

function attach_post_synapse!(self::Neuron, synapse::Segments.Synapse)
    synapticCoordinate = Synapses.get_post_synaptic_coordinate( synapse )
    segmentId, nodeIdInSegment = find_closest_node( self, synapticCoordinate )
    segmentList = get_segment_list(self)
    Segments.attach_post_synapse!(segmentList[segmentId], nodeIdInSegment, synapse) 
end 

"""
    adjust_segment_class(self::Neuron)
adjust the segment class considering the observation that once a segment becomes axon, 
all the children segments should be axon.
"""
function adjust_segment_class!(self::Neuron)
    # adjust class according to the number of pre and post synapses as an initialization 
    for segment in get_segment_list(self)
        Segments.adjust_class!(segment)
    end 
    
    # dendrite parent should all be dendrite 
    for (segmentIndex, segment) in get_segment_list(self) |> enumerate 
        if Segments.DENDRITE_CLASS == segment.class 
            parentSegmentIndex = get_parent_segment_id(self, segmentIndex)
            if parentSegmentIndex > 0
                parentSegment = self[parentSegmentIndex]
                if Segments.AXON_CLASS == parentSegment.class 
                    parentSegment.class = Segments.DENDRITE_CLASS 
                end 
            end 
        end
    end 
    nothing
end

"""
    postprocessing(self::Neuron)
post process after skeletonization
"""
function postprocessing(self::Neuron)
    self = remove_redundent_nodes(self)
    self = remove_hair(self)
    self = remove_subtree_in_soma(self)
    self = remove_terminal_blobs(self)
    self = resample(self, Float32(100))
    self = merge_only_child(self)
    self = smooth(self)
    return self
end  

"""
    remove_segments_with_class(self::Neuron, class::UInt8)
Note that the axon classes should be labeled first!
"""
function remove_segments_with_class(self::Neuron, class::UInt8)
    # we should assume that the axons were labeled to avoid redundent computation 
    #adjust_segment_class!(self)
    segmentIdList = Vector{Int}()
    for (segmentId, segment) in get_segment_list(self) |> enumerate 
        if class == segment.class 
            push!(segmentIdList, segmentId)
        end 
    end 
    remove_segments(self, segmentIdList)
end 

remove_axons(self::Neuron) = remove_segments_with_class(self, Segments.AXON_CLASS)
remove_dendrites(self::Neuron) = remove_segments_with_class(self, Segments.DENDRITE_CLASS)



function downsample_nodes(self::Neuron{T}; nodeNumStep::Int=24) where T 
	@assert nodeNumStep > 1
    segmentList = get_segment_list(self)
    newSegmentList = Vector{Segment{T}}()
    sizehint!(newSegmentList, length(segmentList) )
    for segment in segmentList 
        nodeList = Segments.get_node_list(segment)
        newNodeList = Vector{Segments.Node}()
        for i in 1:nodeNumStep:length(nodeList)
            center = Segments.get_center(segment, 
                                            i: min(i+nodeNumStep-1,length(nodeList)))
            push!(newNodeList, center)
        end 
        newSegment = Segment(newNodeList, class=Segments.get_class(segment))
        push!(newSegmentList, newSegment)
    end
    Neuron{T}(newSegmentList, get_connectivity_matrix(self))
end


"""
    smooth(nodeList::Vector{NTuple{4,T}}, k::Integer) where T 
"""
function smooth(nodeList::Vector{NTuple{4,T}}) where T
    if length(nodeList) < 3 
        return nodeList 
    end  
    newNodeList = Vector{NTuple{4,T}}()
    push!(newNodeList, nodeList[1])
    for i in 2:length(nodeList)-1
        newNode = map((x,y,z)-> (x+y+z) / T(3), nodeList[i-1], nodeList[i], nodeList[i+1])
        push!(newNodeList, newNode)
    end  
    push!(newNodeList, nodeList[end])
    @assert length(newNodeList) == length(nodeList)
    return newNodeList 
end 

"""
    smooth(self::Neuron{T}, k::Integer=3)

make the skeleton smooth. 
use the mean coordinate and radius with two neighboring nodes as new node coordinate 
"""
function smooth(self::Neuron{T}) where T
    newSegmentList = Vector{Segment{T}}()
    segmentList = get_segment_list(self)
    preSynapseList = get_all_pre_synapse_list(self)
    postSynapseList = get_all_post_synapse_list(self)
    
    for (index, segment) in enumerate(segmentList)
        parentSegmentId = get_parent_segment_id(self, index)
        if parentSegmentId > 0 
            parentNode = segmentList[parentSegmentId][end]
            nodeList = [parentNode, Segments.get_node_list(segment) ...]
            newNodeList = smooth(nodeList) 
            # ignore the first node since it is from the parent segment 
            newSegment = Segment(newNodeList[2:end]; class=Segments.get_class(segment))
            push!(newSegmentList, newSegment)
        else 
            # no parent segment
            nodeList = Segments.get_node_list(segment) 
            newNodeList = smooth(nodeList) 
            newSegment = Segment(newNodeList; class=Segments.get_class(segment))
            push!(newSegmentList, newSegment)
        end 
    end 
    @assert length(newSegmentList) == length(segmentList)
    neuron = Neuron{T}(newSegmentList, get_connectivity_matrix(self))

    # reattach synapses 
    attach_pre_synapses!(neuron, preSynapseList)
    attach_post_synapses!(neuron, postSynapseList)
    neuron
end 

"""
    resample( self::Neuron, resampleDistance::Float32 )
resampling the segments by a fixed point distance.
Note that the distance is physical rather than point number.
For example, if the neuron coordinate is based on nm, and the 
resample distance is 1000.0, then the resample distance is 1000 nm.
"""
function resample(self::Neuron{T}, resampleDistance::T) where T 
    newSegmentList = Vector{Segment{T}}()
    segmentList = get_segment_list(self)
    preSynapseList = get_all_pre_synapse_list(self)
    postSynapseList = get_all_post_synapse_list(self)
    for (index, segment) in enumerate(segmentList)
        parentSegmentId = get_parent_segment_id(self, index)
        if parentSegmentId > 0 
            parentNode = segmentList[parentSegmentId][end]
            nodeList = [parentNode, Segments.get_node_list(segment) ...]
            nodeList = resample(nodeList, resampleDistance)
            # ignore the first node since it is from the parent segment 
            newSegment = Segment(nodeList[2:end]; class=Segments.get_class(segment))
            push!(newSegmentList, newSegment)
        else 
            # no parent segment
            nodeList = Segments.get_node_list(segment) 
            nodeList = resample(nodeList, resampleDistance) 
            newSegment = Segment(nodeList; class=Segments.get_class(segment))
            push!(newSegmentList, newSegment)
        end 
    end 
    neuron = Neuron{T}(newSegmentList, get_connectivity_matrix(self))

    # reattach synapses 
    attach_pre_synapses!(neuron, preSynapseList)
    attach_post_synapses!(neuron, postSynapseList)
    neuron
end

function resample(nodeList::Vector{NTuple{4,T}}, resampleDistance::T) where T
    @assert !isempty(nodeList)
    # always keep the first node
    ret = [nodeList[1]]

    # walk through the nodes
    accumulatedNodeDistance = zero(T) 
    @inbounds for index in 1:length(nodeList)-1
        n1 = [nodeList[index  ] ...]
        n2 = [nodeList[index+1] ...]
        nodePairDistance = norm(n1[1:3].-n2[1:3])
        accumulatedNodeDistance += nodePairDistance
        while accumulatedNodeDistance > resampleDistance 
            # walk a step of resample distance 
            # interperate the distance 
            interpretRatio = 1 - (accumulatedNodeDistance-resampleDistance) / nodePairDistance
            @assert (interpretRatio >= zero(T) && interpretRatio <= one(T)) || 
                    isapprox(interpretRatio, 0; atol=6)
            interpretedNode = n1 .+ (n2.-n1) .* interpretRatio 
            push!(ret, (interpretedNode...,))

            # update accumulated node distance  
            accumulatedNodeDistance -= resampleDistance 
        end 
    end
    if accumulatedNodeDistance > zero(T)
        # if there is remaining distance, add the last node
        push!(ret, nodeList[end])
    end
    @assert !isempty(ret)
    ret 
end

"""
    build_tree_with_position(self::Neuron{T}; leafSize=1) where T 
    
Return: 
    kdTree::KDTree, the built kdtree for node query 
    nodePosition::Matrix{Int}, the first row is segment id, 
            the second row is node id in the segment.
            the column number is the node number of this neuron.
"""
function build_tree_with_position(self::Neuron{T}; leafSize::Integer=1) where T 
    N = get_node_num(self)

    nodeCoord = Matrix{T}(undef, 3, N)
    # the node position in the neuron
    # the first row is segment id, and the second row is the node id in segment 
    nodePosition = Matrix{Int}(undef, 2, N)

    k = 0
    for (segmentId,segment) in enumerate(self)
        for (nodeId, node) in enumerate(segment)
            k += 1
            nodeCoord[:,k] = [node[1:3]...]
            nodePosition[1,k] = segmentId 
            nodePosition[2,k] = nodeId 
        end
    end 
    kdTree = KDTree(nodeCoord; leafsize=leafSize)
    return kdTree, nodePosition 
end

function find_closest_node(self::Neuron{T}, seedNode::NTuple{3,T}) where T
    kdTree, nodePosition = build_tree_with_position(self)
    find_closest_node(kdTree, nodePosition, seedNode)
end 

function find_closest_node(kdTree::KDTree, nodePosition::Matrix{Int}, 
                           seedNode::Union{Tuple, Vector})
    idxs, dists = knn(kdTree, [seedNode[1:3]...], 1, false)
    segmentId, nodeIdInSegment = nodePosition[:, idxs[1]]
    return segmentId, nodeIdInSegment, dists[1]
end 

"""
    find_closest_node_V1(self::Neuron{T}, seedNode::NTuple{3,T}) where T

Return:
    closestSegmentId::Int, the closest segment Id 
    closestNodeIdInSegment::Int, the closest node id in the segment  
    distance::T, the closest distance 
"""
function find_closest_node_V1(self::Neuron{T}, seedNode::NTuple{3,T}) where T 
    segmentList = get_segment_list(self)
    
    # initialization 
    closestSegmentId = 0
    closestNodeIdInSegment = 0
    distance = typemax(T)
    
    @assert !isempty(segmentList)
    
    for (segmentId, segment) in enumerate(segmentList)
        # distance from bounding box give the maximum bound of closest distance
        # this was used to filter out far away segments quickly
        # no need to compute all the distances between nodes
        # this is the lower bound of the distance
        bbox_distance = Segments.get_bounding_box_distance(segment, seedNode)
        if bbox_distance < distance 
            d, nodeIdInSegment = Segments.distance_from(segment, seedNode )
            if d < distance 
                distance = d
                closestSegmentId = segmentId 
                closestNodeIdInSegment = nodeIdInSegment 
            end 
        end 
    end
    @assert closestNodeIdInSegment > 0 
    @assert closestNodeIdInSegment <= length(segmentList[closestSegmentId])
    closestSegmentId, closestNodeIdInSegment, distance 
end 


"""
    find_merging_terminal_node_id(self::Neuron, seedNode2::NTuple{4, Float32}) 

# Note that the weighted version was commented out since it seems not working well
#the distance was weighted according to angle θ: d' = d * tan(θ)
#this function has a nice property that shrinking the distance within 45 degrees  
#and enlarge it with angle greater than 45 degrees 
"""
function find_merging_terminal_node_id(self::Neuron{T}, seedNode2::NTuple{4, T}) where T
    # initialization 
    closestTerminalSegmentId1 = 0
    terminalNodeIdInSegment1 = 0
    weightedDistance = typemax(T)

    terminalSegmentIdList1 = get_terminal_segment_id_list(self)
    segmentList1 = get_segment_list(self)

    for terminalSegmentId1 in terminalSegmentIdList1 
        terminalSegment1 = segmentList1[terminalSegmentId1]
        terminalNode1 = terminalSegment1[end]
        b = [ [seedNode2[1:3]...] .- [terminalNode1[1:3]...] ]
        physicalDistance::T = norm(b)
        
        if physicalDistance < weightedDistance 
            weightedDistance = physicalDistance 
            closestTerminalSegmentId1 = terminalSegmentId1
            terminalNodeIdInSegment1 = length(terminalSegment1)
        end 
        
       # if length(terminalSegment1) > 1
       #     terminalNode0 = terminalSegment1[end-1]
       # else 
       #     parentSegmentId = get_parent_segment_id(self, terminalSegmentId1)
       #     terminalNode0 = segmentList1[parentSegmentId][end]
       # end 

       # # get angle θ
       # a = [ [terminalNode1[1:3]...] .- [terminalNode0[1:3]...] ]
       # # sometimes the value could be larger than 1.0 or less than -1.0 due to numerical precision stabilities
       # dp = dot(a,b) / (norm(a)*norm(b))
       #         
       # # euclidean distance
       # #if theta > pi/2
       # if dp < zero(Float32)
       #     # if we use one(Float32)), the tan value will become positive due to numerical precision. 
       #     # acos(-one(Float32)) |> tan = 8.7422784f-8 
       #     dp = max(-0.9999999f0, dp)
       #     theta = acos(dp)
       #     #@show theta

       #     distance = physicalDistance * (-tan(theta))
       #     #@assert distance > zero(typeof(distance))
       #     if distance < weightedDistance
       #         weightedDistance = distance 
       #         closestTerminalSegmentId1 = terminalSegmentId1
       #     end
       # end 
    end

    return closestTerminalSegmentId1, terminalNodeIdInSegment1, weightedDistance  
end 

"""
    find_seed_node_id(neuron::Neuron, nodeNet::NodeNets, collectedFlagVec::BitArray{1})

find the closest terminal node in uncollected node set as new growing seed.  
"""
function find_seed_node_id(neuron::Neuron, nodeNet::NodeNet, collectedFlagVec::BitArray{1})
    # number 1 means the alread grown neuron 
    # number 2 means the id in the raw node net 
    segmentList1 = get_segment_list(neuron)
    # the new seed will be chosen this terminal node set 
    terminalNodeIdList2 = NodeNets.get_terminal_node_id_list(nodeNet)

    # initialization
    seedNodeId2 = 0
    distance = typemax(Float32)
    
    @assert !isempty(segmentList1)
    @assert !all(collectedFlagVec)
    
    for candidateSeedId2 in terminalNodeIdList2 
        if !collectedFlagVec[candidateSeedId2]
            # only choose from the uncollected ones 
            node2 = nodeNet[candidateSeedId2] 
            for (segmentId1, segment1) in enumerate(segmentList1)
                # distance from bounding box give the maximum bound of closest distance
                # this was used to filter out far away segments quickly
                # no need to compute all the distances between nodes
                # this is the lower bound of the distance
                bbox_distance = Segments.get_bounding_box_distance(segment1, node2)
                if bbox_distance < distance 
                    d, _ = Segments.distance_from(segment1, node2)
                    if d < distance 
                        distance = d
                        seedNodeId2 = candidateSeedId2 
                    end 
                end 
            end
        end 
    end
    return seedNodeId2 
end 


"""
    find_closest_node(neuron::Neuron, nodeNet::NodeNet, collectedFlagVec::BitArray{1})

find the uncollected node which is closest to the segment list 
"""
function find_closest_node(neuron::Neuron, nodeNet::NodeNet, collectedFlagVec::BitArray{1})
    segmentList = get_segment_list(neuron)
    nodes = NodeNets.get_node_list(nodeNet)

    # initialization 
    closestNodeId = 0
    mergingSegmentId = 0
    mergingNodeIdInSegment = 0
    distance = typemax(Float32)
    
    @assert !isempty(segmentList)
    @assert length(nodes) == length(collectedFlagVec)
    @assert !all(collectedFlagVec)
    
    for (nodeId, node) in enumerate(nodes)
        if !collectedFlagVec[nodeId] && NodeNets.is_terminal_node(nodeNet, nodeId)
            for (segmentId, segment) in enumerate(segmentList)
                # distance from bounding box give the maximum bound of closest distance
                # this was used to filter out far away segments quickly
                # no need to compute all the distances between nodes
                # this is the lower bound of the distance
                bbox_distance = Segments.get_bounding_box_distance(segment, node)
                if bbox_distance < distance 
                    d, nodeIdInSegment = Segments.distance_from(segment, node)
                    if d < distance 
                        distance = d
                        closestNodeId = nodeId 
                        mergingSegmentId = segmentId 
                        mergingNodeIdInSegment = nodeIdInSegment 
                    end 
                end 
            end
        end 
    end
    # the segmentId should be inside the segmentList
    #@assert mergingSegmentId > 0
    #@assert mergingSegmentId <= length(segmentList)
    return closestNodeId, mergingSegmentId, mergingNodeIdInSegment 
end 

function add_offset(self::Neuron{T}, offset::Union{Vector, Tuple}) where T
    segmentList = Vector{Segment{T}}()
    for segment in self.segmentList 
        newSegment = Segments.add_offset(segment, offset)
        push!(segmentList, newSegment)
    end 
    Neuron{T}(segmentList, self.connectivityMatrix)
end 

end # module
