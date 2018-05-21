module Neurons
include("Segments.jl")
using .Segments
using ..NodeNets
using ..SWCs
using ..Utils.BoundingBoxes 
using LsqFit
using ImageFiltering
using OffsetArrays
using DataFrames
using .Segments.Synapses 
using Alembic


const ONE_UINT32 = UInt32(1)
const ZERO_FLOAT32 = Float32(0)
const EXPANSION = (ONE_UINT32, ONE_UINT32, ONE_UINT32)
const VOXEL_SIZE = (1000, 1000, 1000)

const DOWNSAMPLE_NODE_NUM_STEP = 24

export Neuron

mutable struct Neuron 
    # x,y,z,r, the coordinates should be in physical coordinates
    segmentList ::Vector{Segment}
    connectivityMatrix ::SparseMatrixCSC{Bool, Int}
end 

"""
    Neuron
a neuron modeled by interconnected segments 
"""
function Neuron(nodeNet::NodeNet)
    # the properties from nodeNet
    nodes = NodeNets.get_node_list(nodeNet)
    radii = NodeNets.get_radii(nodeNet)
    # flags labeling whether this node was collected to the net
    collectedFlagVec = zeros(Bool, length(nodes))
    # connectivity matrix of nodes in nodeNet 
    nodesConnectivityMatrix  = NodeNets.get_connectivity_matrix(nodeNet)
    # locate the root node with largest radius
    # theoritically this should be the center of soma 
    _, rootNodeIndex = findmax(radii)
    seedNodeIndexList::Vector = [rootNodeIndex]
    # grow the main net
    neuron = Neuron!(rootNodeIndex, nodeNet, collectedFlagVec)
    #return neuron # only use the main segment for evaluation

    while !all(collectedFlagVec)
        # there exist some uncollected nodes 
        # find the uncollected node that is closest to the terminal nodes as seed
        # terminal node should have high priority than normal nodes. 
        seedNodeIndex, closestSegmentIndex, closestDistance = 
            find_closest_terminal_node_index(neuron, nodeNet, collectedFlagVec)
        if closestDistance > 6000
            # for some spine like structures, they should connect to nearest node
            # rather than terminal point 
            seedNodeIndex, closestSegmentIndex = 
                find_closest_node_index(neuron, nodeNet, collectedFlagVec)
        end 

        # get the terminal node index of the terminal segment
        segmentList = get_segment_list(neuron)
        @assert closestSegmentIndex <= length(segmentList)
        closestSegment = segmentList[closestSegmentIndex]
        # the last node should be the terminal node
        closestNodeIndexInSegment = length(closestSegment)

        # grow a new net from seed, and mark the collected nodes 
        subnet = Neuron!(seedNodeIndex, nodeNet, collectedFlagVec) 
        # merge the subnet to main net 
        neuron = merge( neuron, subnet, 
                        closestSegmentIndex, closestNodeIndexInSegment)
    end 
    neuron     
end 

"""
    Neuron!(seedNodeIndex::Integer, nodeNet::NodeNet, collectedFlagVec::Vector{Bool})

build a net from a seed node using connected component
mark the nodes in this new net as collected, so the collectedFlagVec was changed.
"""
function Neuron!(seedNodeIndex::Integer, nodeNet::NodeNet, 
                                            collectedFlagVec::Vector{Bool})
    # initialization
    segmentList = Vector{Segment}()
    parentSegmentIndexList = Vector{Int}()
    childSegmentIndexList  = Vector{Int}()

    nodes = NodeNets.get_node_list(nodeNet)
    nodesConnectivityMatrix = NodeNets.get_connectivity_matrix( nodeNet )
    
    # the seed of segment should record both seed node index 
    # and the the parent segment index
    # the parent segment index of root segment is -1 
    segmentSeedList = [(seedNodeIndex, -1)]
    # depth first search
    while !isempty(segmentSeedList)
        seedNodeIndex, segmentParentIndex = pop!(segmentSeedList)
        # grow a segment from seed
        nodeListInSegment = Vector{NTuple{4, Float32}}()
        seedNodeIndexList = [seedNodeIndex]
        while true
            # construct this segment
            seedNodeIndex = pop!(seedNodeIndexList)
            # push the seed node
            push!(nodeListInSegment, nodes[seedNodeIndex])
            # label this node index as collected
            collectedFlagVec[ seedNodeIndex ] = true 

            # find the connected nodes
            connectedNodeIndexList,_ = findnz(nodesConnectivityMatrix[:, seedNodeIndex])
            # exclude the collected nodes
            connectedNodeIndexList = connectedNodeIndexList[ .!collectedFlagVec[connectedNodeIndexList] ] 

            if length(connectedNodeIndexList) == 1
                # belong to the same segment
                push!(nodeListInSegment, nodes[ connectedNodeIndexList[1] ])
                push!(seedNodeIndexList, connectedNodeIndexList[1])
            else
                # terminal branching point or multiple branching points
                # finish constructing this segment
                segment = Segment(nodeListInSegment)
                push!(segmentList, segment)
                if segmentParentIndex != -1
                    # this is not the root segment, establish the segment connection
                    push!(parentSegmentIndexList, segmentParentIndex)
                    push!(childSegmentIndexList,  length(segmentList))
                end 
                # seed new segments
                # if this is terminal segment, no seed will be pushed
                for index in connectedNodeIndexList 
                    push!(segmentSeedList, (index, length(segmentList)))
                end 
                break
            end 
        end 
    end
    @assert length(parentSegmentIndexList) == length(childSegmentIndexList)
    # note that the connectivity matrix should be a square matrix for easier use
    connectivityMatrix = spzeros(Bool, length(segmentList), length(segmentList))
    tempConnectivityMatrix = sparse(parentSegmentIndexList, childSegmentIndexList, 
                                ones(Bool,length(childSegmentIndexList)))
    connectivityMatrix[1:size(tempConnectivityMatrix,1), 
                       1:size(tempConnectivityMatrix,2)] = tempConnectivityMatrix
    for segment in segmentList 
        @assert !isempty(segment)
    end
    Neuron(segmentList, connectivityMatrix)
end 

function Neuron( seg::Array{T,3}; obj_id::T = convert(T,1), 
                     expansion::NTuple{3,UInt32}=EXPANSION ) where T
    nodeNet = NodeNet( seg; obj_id = obj_id, expansion = expansion )
    Neuron( nodeNet )
end

function Neuron( swc::SWC )
    nodeNet = NodeNet( swc )
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
function load_swc( fileName::AbstractString )
    swcString = readstring( fileName )
    Neuron( swcString )
end

function load_swc_bin( fileName::AbstractString )
    swc = SWCs.load_swc_bin( fileName )
    Neuron( swc )
end 

function save(self::Neuron, fileName::AbstractString)
    swc = SWC(self)
    SWCs.save( swc, fileName )
end 
save_swc = save

function save_swc_bin( self::Neuron, fileName::AbstractString )
    SWCs.save_swc_bin( SWCs.SWC(self), fileName )
end 

####################### properties ##############

"""
    get_root_segment_index( self::Neuron )

the first segment should be the root segment, and the first node should be the root node 
"""
function get_root_segment_index(self::Neuron) 1 end 
function get_root_segment(self::Neuron) 
    get_segment_list(self)[get_root_segment_index(self)]
end 

"""
    get_root_node( self::Neuron )
the returned root node is a tuple of Float64, which represents x,y,z,r
"""
function get_root_node( self::Neuron ) 
    rootSegmentIndex = get_root_segment_index(self)
    rootSegment = get_segment_list(self)[ rootSegmentIndex ]
    rootSegment[1]
end 

"""
    get_num_segments(self::Neuron)
"""
function get_num_segments(self::Neuron) length(self.segmentList) end

"""
    get_num_branching_points(self::Neuron)
"""
function get_num_branching_points(self::Neuron)
    numSegmentingPoint = 0
    for index in 1:get_num_segments(self)
        childrenSegmentIndexList = get_children_segment_index_list(self, index)
        if length(childrenSegmentIndexList) > 0
            numSegmentingPoint += 1
        end 
    end 
    numSegmentingPoint
end 

function get_segment_list(self::Neuron) self.segmentList end 
function get_connectivity_matrix(self::Neuron) self.connectivityMatrix end 

"""
    get_segment_order_list( self::Neuron )
following the Julia indexing style, the root segment order is 1.
"""
function get_segment_order_list(self::Neuron)
    segmentOrderList = Vector{Int}()
    sizehint!(segmentOrderList, get_num_segments(self))
    
    index2order = Dict{Int,Int}()
    # the first one is segment index, the second one is segment order 
    seedSegmentIndexOrderList = Vector{NTuple{2,Int}}([(1,1)])
    while !isempty( seedSegmentIndexOrderList )
        segmentIndex, segmentOrder = pop!(seedSegmentIndexOrderList)
        index2order[segmentIndex] = segmentOrder 
        childrenSegmentIndexList = get_children_segment_index_list(self, segmentIndex)
        for childSegmentIndex in childrenSegmentIndexList 
            push!(seedSegmentIndexOrderList, (childSegmentIndex, segmentOrder+1))
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
function get_segment_length_list( self::Neuron ) map(length, get_segment_list(self)) end

"""
    get_node_num( self::Neuron )
get total number of nodes  
"""
function get_node_num( self::Neuron )
    sum(map(length, get_segment_list(self)))
end 

"""
    get_node_list(self::Neuron)
get the node list. the first one is root node.
"""
function get_node_list( self::Neuron )
    nodeList = Vector{NTuple{4, Float32}}()
    sizehint!(nodeList, get_node_num(self))
    for segment in get_segment_list(self)
        append!(nodeList, Segments.get_node_list(segment))
    end 
    nodeList
end 

"""
    get_edge_list( self::Neuron )
get the edges with type of Vector{NTuple{2, Int}}
"""
function get_edge_list( self::Neuron )
    edgeList = Vector{NTuple{2,Int}}()
    segmentStartNodeIndexList = Vector{Int}()
    segmentStopNodeIndexList = Vector{Int}()
    # total number of nodes
    nodeNum = 0
    for (segmentIndex, segment) in enumerate(get_segment_list(self))
        push!(segmentStartNodeIndexList, nodeNum+1)
        # build the edges inside segment
        for nodeIndex in nodeNum+1:nodeNum+length(segment)-1
            push!(edgeList, (nodeIndex, nodeIndex+1))
        end 
        # update node number
        nodeNum += length(segment)
        push!(segmentStopNodeIndexList,  nodeNum)
    end 
    # add segment connections
    parentSegmentIndexList, childSegmentIndexList, _ = findnz( get_connectivity_matrix(self) )
    for (index, parentSegmentIndex) in enumerate( parentSegmentIndexList )
        childSegmentIndex = childSegmentIndexList[ index ]
        parentNodeIndex = segmentStopNodeIndexList[ parentSegmentIndex ]
        childNodeIndex  = segmentStartNodeIndexList[ childSegmentIndex ]
        push!( edgeList, (parentNodeIndex, childNodeIndex) )
    end 
    edgeList 
end 

function get_path_to_soma_length(self::Neuron, synapse::Synapse)
    closestSegmentIndex, closestNodeIndex = find_closest_node_index( self, synapse )
    get_path_to_soma_length( self, closestSegmentIndex; nodeIndex=closestNodeIndex )
end 

function get_path_to_soma_length(self::Neuron, segmentIndex::Integer; 
                        segmentList::Vector{Segment}=get_segment_list(self), 
                        nodeIndex::Int = length(segmentList[segmentIndex]), 
                        segmentPathLengthList::Vector = get_segment_path_length_list(self))
    path2RootLength = Segments.get_path_length( segmentList[segmentIndex]; 
                                                        nodeIndex=nodeIndex )
    while true 
        parentSegmentIndex = get_parent_segment_index(self, segmentIndex )
        if parentSegmentIndex < 1 
            # root segment do not have parent 
            break 
        else
            path2RootLength += segmentPathLengthList[ parentSegmentIndex ]
            segmentIndex = parentSegmentIndex 
        end 
    end
    path2RootLength 
end

function get_pre_synapse_to_soma_path_length_list(self::Neuron; 
                        segmentList::Vector=get_segment_list(self),
                        segmentPathLengthList::Vector=get_segment_path_length_list(self))
    preSynapseToSomaPathLengthList = Vector{Float32}()
    for (segmentIndex, segment) in enumerate( segmentList )
        preSynapseList = Segments.get_pre_synapse_list( segment )
        nodeIndexList, _ = findnz( preSynapseList )
        pathToSomaLengthList = map(nodeIndex->get_path_to_soma_length(self, segmentIndex; 
                                    segmentList=segmentList,
                                    segmentPathLengthList=segmentPathLengthList,
                                    nodeIndex=nodeIndex), nodeIndexList)
        append!(preSynapseToSomaPathLengthList, pathToSomaLengthList)
    end
    preSynapseToSomaPathLengthList 
end 

function get_post_synapse_to_soma_path_length_list(self::Neuron; 
                        segmentList::Vector=get_segment_list(self),
                        segmentPathLengthList::Vector=get_segment_path_length_list(self))
    postSynapseToSomaPathLengthList = Vector{Float32}()
    for (segmentIndex, segment) in enumerate( segmentList )
        postSynapseList = Segments.get_post_synapse_list( segment )
        nodeIndexList, _ = findnz( postSynapseList )
        pathToSomaLengthList = map(nodeIndex->get_path_to_soma_length(self, segmentIndex; 
                                    segmentList=segmentList,
                                    segmentPathLengthList=segmentPathLengthList,
                                    nodeIndex=nodeIndex), nodeIndexList)
        append!(postSynapseToSomaPathLengthList, pathToSomaLengthList)
    end
    postSynapseToSomaPathLengthList 
end 


"""
    get_segment_path_length_list(self::Neuron)
get euclidean path length of each segment 
"""
function get_segment_path_length_list( self::Neuron; segmentList = get_segment_list(self) )
    ret = Vector{Float64}()
    for (index, segment) in enumerate( segmentList )
        segmentPathLength = Segments.get_path_length( segment )
        # add the edge length to parent node
        parentSegmentIndex = get_parent_segment_index(self, index)
        if parentSegmentIndex > 0
            parentSegment = segmentList[ parentSegmentIndex ]
            parentNode = parentSegment[ end ]
            node = segment[1]
            segmentPathLength += norm( [map((x,y)->x-y, node[1:3], parentNode[1:3])...] )
        end 
        push!(ret, segmentPathLength)
    end 
    ret 
end 

function get_node_distance_list(neuron::Neuron)
    nodeList = Neurons.get_node_list(neuron)
    edgeList = Neurons.get_edge_list(neuron)
    nodeDistanceList = Vector{Float32}()
    sizehint!(nodeDistanceList, length(edgeList))
    for (src, dst) in edgeList
        d = norm( [map(-, nodeList[src][1:3], nodeList[dst][1:3])...] )
        push!(nodeDistanceList, d)
    end 
    nodeDistanceList
end 

function get_total_path_length(self::Neuron)
    get_segment_path_length_list(self) |> sum
end 

"""
    get_segment_end_node_index_list( self::Neuron )

get a vector of integer, which represent the node index of segment end 
"""
function get_segment_end_node_index_list( self::Neuron ) cumsum( get_segment_length_list(self) ) end 

"""
    get_num_nodes(self::Neuron)

get number of nodes 
"""
function get_num_nodes( self::Neuron )
    segmentList = get_segment_list(self)
    sum(map(length, segmentList))
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
    get_children_segment_index_list(self::Neuron, parentSegmentIndex::Integer)
"""
function get_children_segment_index_list(self::Neuron, parentSegmentIndex::Integer)
    childrenSegmentIndexList,_ = findnz(self.connectivityMatrix[parentSegmentIndex, :])
    childrenSegmentIndexList 
end

function get_parent_segment_index( self::Neuron, childSegmentIndex::Integer )
    parentSegmentIndexList,_ = findnz(self.connectivityMatrix[:, childSegmentIndex])
    @assert length(parentSegmentIndexList) <= 1
    if isempty( parentSegmentIndexList ) 
        # no parent, this is a root segment
        return 0
    else 
        return parentSegmentIndexList[1]
    end 
end

"""
    get_subtree_segment_index_list( self, segmentInde )
get the segment index list of subtree 
"""
function get_subtree_segment_index_list( self::Neuron, segmentIndex::Integer )
    @assert segmentIndex > 0 && segmentIndex <= get_num_segments(self)
    subtreeSegmentIndexList = Vector{Integer}()
    seedSegmentIndexList = [segmentIndex]
    while !isempty( seedSegmentIndexList )
        seedSegmentIndex = pop!( seedSegmentIndexList )
        push!(subtreeSegmentIndexList, seedSegmentIndex)
        childrenSegmentIndexList = get_children_segment_index_list( self, seedSegmentIndex )
        append!(seedSegmentIndexList, childrenSegmentIndexList)
    end 
    subtreeSegmentIndexList 
end 

function get_terminal_segment_index_list( self::Neuron; startSegmentIndex::Integer = 1 )
    terminalSegmentIndexList = Vector{Int}()
    segmentIndexList = [ startSegmentIndex ]
    for segmentIndex in segmentIndexList 
        childrenSegmentIndexList = get_children_segment_index_list( self, segmentIndex )
        if isempty( childrenSegmentIndexList ) 
            # this is a terminal segment
            push!(terminalSegmentIndexList, segmentIndex )
        else
            append!(segmentIndexList, childrenSegmentIndexList)
        end 
    end 
    terminalSegmentIndexList 
end

function get_terminal_node_list( self::Neuron; startSegmentIndex::Integer = 1 )
    terminalSegmentIndexList = get_terminal_segment_index_list( self )
    map( x -> Segments.get_node_list(x)[end], get_segment_list(self)[ terminalSegmentIndexList ] )
end 

"""
    get_branching_angle( self::Neuron, segmentIndex::Integer; nodeDistance::AbstractFloat  )
if the node is too close the angle might not be accurate. For example, if they are voxel neighbors, the angle will alwasys be 90 or 45 degree. Thus, we have nodeDistance here to help make the nodes farther away.
Note that the acos returns angle in the format of radiens.
"""
function get_branching_angle( self::Neuron, segmentIndex::Integer; nodeDistance::Real = 5000.0 )
    segment = self[segmentIndex]
    parentSegmentIndex = get_parent_segment_index(self, segmentIndex)
    if parentSegmentIndex < 1
        # this is root segment
        return 0.0
    end
    parentSegment = self[ get_parent_segment_index(self, segmentIndex) ]
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
        warn("parent node is the same with branching node: $(branchingNode)")
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
        warn("child node is the same with branching node: $(branchingNode)")
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
    for (centerIndex, center) in enumerate(diskCenterList)
        for (radiusIndex, radius) in enumerate(radiusList)
            for node in nodeList 
                if Segments.get_nodes_distance( center,  node) < radius
                    averageMassList[radiusIndex] += 1.0 / length(diskCenterList)
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
        push!(voxelSet, (voxelCoordinate...))
    end 
    boundingBox = Segments.BoundingBox( voxelSet )
    range = BoundingBoxes.get_unit_range(boundingBox)
    sz = size(boundingBox)
    # initialize the map 
    mask = zeros(Bool, sz)
    mask = OffsetArray(mask, range...)

    @unsafe for voxel in voxelSet
        # add a small number to make sure that there is no 0
        mask[voxel...] = true 
    end 
    mask
end 

function get_2d_binary_projection(self::Neuron; axis=3, voxelSize=VOXEL_SIZE)
    nodeList = get_node_list(self)
    voxelSet = Set{NTuple{3, Int}}()
    for node in nodeList 
        voxelCoordinate = map((x,y)->round(Int, x/y), node[1:3], voxelSize)
        push!(voxelSet, (voxelCoordinate...))
    end 
    boundingBox = Segments.BoundingBox( voxelSet ) 
    range = Segments.BoundingBoxes.get_unit_range(boundingBox)[1:2]
    sz = size(boundingBox)[1:2]
    # initialize the map                                                   
    mask = zeros(Bool, sz)
    mask = OffsetArray(mask, range...)

    @unsafe for voxel in voxelSet 
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
    kernel = fill(gaussianFilterStd, 3) |> Kernel.gaussian
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
    newRange = map((x,y,s)->(x-round(Int,y/s)), 
                   indices(densityMap), somaCorrdinate, voxelSize)
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
    intersectIndices = map(intersect, indices(self), indices(other))
    d = 0.0
    # add the distance of intersection
    d += sum(abs2.(self[intersectIndices...] .- other[intersectIndices...])) 
    # add the distance of self from non-intersect regions
    selfDiff  = map(setdiff, indices(self),  intersectIndices)
    otherDiff = map(setdiff, indices(other), intersectIndices)
    d += sum(abs2.(self[selfDiff...]))
    d += sum(abs2.(other[otherDiff...]))
    d
end 


function _cross_correlation_along_axis(self::Array{T,3}, other::Array{T,3}, 
                                                                axis::Int) where {T}
    projection1 = maximum(self,  axis)
    projection2 = maximum(other, axis)
	projection1 = squeeze(projection1, axis)
	projection2 = squeeze(projection2, axis)
    
    
    template = projection2 
    image = projection1  
    if length(template) > length(image)
        image, template = template, image 
    end
    
    # pad the image with mean value to make sure that the template is smaller than image
    if size(image, 1) < size(template, 1)
        padSize = size(template,1) - size(image,1)
        image = ImageFiltering.padarray( image, Fill(mean(image), (padSize,0), (0,0)) )
    elseif size(image,2) < size(template,2)
        padSize = size(template,2) - size(image,2)
        image = ImageFiltering.padarray( image, Fill(mean(image), (0,padSize), (0,0)))
    end 
    
    # projection1 is template, and projection2 is the image
    return Alembic.normxcorr2_preallocated(template|>parent, image|>parent; shape="same")
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
function Base.getindex(self::Neuron, index::Integer)
    get_segment_list(self)[index]
end

"""
merge two nets at a specific segment location.
the root node of second net will be connected to the first net
"""
function Base.merge(self::Neuron, other::Neuron, 
                    closestSegmentIndex::Integer, closestNodeIndexInSegment::Integer)
    @assert !isempty(self)
    @assert !isempty(other)
    @assert closestSegmentIndex > 0
    segmentList1 = get_segment_list( self  )
    segmentList2 = get_segment_list( other )
    num_segments1 = get_num_segments(self)
    num_segments2 = get_num_segments(other)
    
    if closestNodeIndexInSegment == length(segmentList1[closestSegmentIndex])
        println("connecting to a segment end...")
        childrenSegmentIndexList = get_children_segment_index_list(self, closestSegmentIndex)
        if length(childrenSegmentIndexList) > 0
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
            mergedConnectivityMatrix[closestSegmentIndex, num_segments1+1] = true
            return Neuron(mergedSegmentList, mergedConnectivityMatrix)
        else 
            println("the closest segment is a terminal segment, merge the root segment of the other net")
            if num_segments2 == 1
                #  only one segment of the other net
                mergedSegmentList = copy(segmentList1)
                mergedSegmentList[closestSegmentIndex] = 
                            merge(mergedSegmentList[closestSegmentIndex], segmentList2[1])
                return Neuron(mergedSegmentList, self.connectivityMatrix)
            else 
                # the connection point is the end of a segment, no need to break segment
                # merge the root segment of second net to the first net
                # assume that the first segment is the root segment
                mergedSegmentList = vcat( segmentList1, segmentList2[2:end] )
                mergedSegmentList[closestSegmentIndex] = 
                            merge(mergedSegmentList[closestSegmentIndex], segmentList2[1])

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
                childrenSegmentIndexList2 = get_children_segment_index_list(other, 1)
                for childSegmentIndex2 in childrenSegmentIndexList2
                    mergedConnectivityMatrix[ closestSegmentIndex, 
                                              num_segments1+childSegmentIndex2-1 ] = true 
                end
                return Neuron(mergedSegmentList, mergedConnectivityMatrix)
            end 
        end 
    else 
        #println("need to break the closest segment and rebuild connectivity matrix")
        total_num_segments = num_segments1 + 1 + num_segments2 
        mergedSegmentList = segmentList1 
        mergedConnectivityMatrix = spzeros(Bool, total_num_segments, total_num_segments)
        
        # need to break the segment and then stitch the new segments
        segmentPart1, segmentPart2 = split(segmentList1[closestSegmentIndex], 
                                                    closestNodeIndexInSegment)
        mergedSegmentList[closestSegmentIndex] = segmentPart1 
        mergedConnectivityMatrix[1:size(self.connectivityMatrix,1), 
                                 1:size(self.connectivityMatrix,2)] = 
                                                        self.connectivityMatrix 
        # reconnect the breaked two segments
        push!(mergedSegmentList, segmentPart2)
        mergedConnectivityMatrix[closestSegmentIndex, num_segments1+1] = true 

        # redirect the children segments to segmentPart2
        childrenSegmentIndexList = get_children_segment_index_list(self, closestSegmentIndex)
        for childSegmentIndex in childrenSegmentIndexList
            # remove old connection
            mergedConnectivityMatrix[closestSegmentIndex, childSegmentIndex] = false
            # build new connection
            mergedConnectivityMatrix[num_segments1+1, childSegmentIndex] = true 
        end 

        # merge the other net
        mergedSegmentList = vcat(mergedSegmentList, segmentList2)
        mergedConnectivityMatrix[
                num_segments1+2 : num_segments1+2+size(other.connectivityMatrix,1)-1, 
                num_segments1+2 : num_segments1+2+size(other.connectivityMatrix,2)-1] = 
                                                                other.connectivityMatrix

        # establish the connection between two nets
        mergedConnectivityMatrix[closestSegmentIndex, num_segments1+2] = true

        # create new merged net
        return Neuron(mergedSegmentList, mergedConnectivityMatrix)
    end
end 

function Base.isempty(self::Neuron)    isempty(self.segmentList) end 

"""
    Base.split(self::Neuron, splitSegmentIndex::Integer; nodeIndexInSegment::Integer=0)

split the segment net into two neurons
the nodeIndexInSegment will be included in the first main net including the original root node  
"""
function Base.split(self::Neuron, splitSegmentIndex::Integer; nodeIndexInSegment::Integer=1)
    if splitSegmentIndex == 1 && nodeIndexInSegment < 1 
        return Neuron(), self
    end 

    segmentList = get_segment_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    segmentIndexListInSubtree2 = get_subtree_segment_index_list( self, splitSegmentIndex )
    subtree1Segment, subtree2RootSegment = split(segmentList[splitSegmentIndex], 
                                                                nodeIndexInSegment)

    # build the split out subtree2, which do not contain the root node
    segmentList2 = [subtree2RootSegment]
    append!(segmentList2, segmentList[ segmentIndexListInSubtree2 ])
    segmentIndexMap = Dict{Int,Int}()
    sizehint!(segmentIndexMap, length(segmentIndexListInSubtree2))
    for (segmentIndexInSubtree2, mainTreeSegmentIndex ) in 
                                                    enumerate(segmentIndexListInSubtree2)
        segmentIndexMap[ mainTreeSegmentIndex ] = segmentIndexInSubtree2 
    end
    connectivityMatrix2 = spzeros( Bool, length(segmentList2), length(segmentList2) )
    for (segmentIndexInSubtree2, segmentIndexInOriginalTree) in 
                                            enumerate(segmentIndexListInSubtree2)
        parentSegmentIndexInOriginalTree = 
                                get_parent_segment_index(self, segmentIndexInOriginalTree)
        if haskey(segmentIndexMap, parentSegmentIndexInOriginalTree)
            parentSegmentIndexInSubtree2 = segmentIndexMap[ parentSegmentIndexInOriginalTree ]
            connectivityMatrix2[ parentSegmentIndexInSubtree2, segmentIndexInSubtree2 ] = true 
        else 
            # this is the root segment of new tree 
            println("find root of new tree in the original tree : $(parentSegmentIndexInOriginalTree)")
        end 
    end 
    subtree2 = Neuron( segmentList2, connectivityMatrix2 )

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    segmentList1 = Vector{Segment}()
    segmentIndexListInSubtree1 = Vector{Int}()
    segmentIndexMap = Dict{Int,Int}()
    # the first one is the root segment
    seedSegmentIndexList = [1]
    # build segmentList1 without the splitout segment, will deal with it later
    while !isempty( seedSegmentIndexList )
        seedSegmentIndex = pop!( seedSegmentIndexList )
        if seedSegmentIndex != splitSegmentIndex 
            push!(segmentList1, segmentList[seedSegmentIndex])
            push!(segmentIndexListInSubtree1, seedSegmentIndex)
            # map the old segment index to the new segment index
            # 3 => 2 means the 3rd segment of old tree is the 2nd segment of new tree
            segmentIndexMap[ seedSegmentIndex ] = length(segmentList1)
            childrenSegmentIndexList = get_children_segment_index_list(self, seedSegmentIndex)
            append!(seedSegmentIndexList, childrenSegmentIndexList)
        else 
            # deal with the split out segment
            if nodeIndexInSegment != 1
                @assert length(subtree1Segment) != 0
                push!(segmentList1, subtree1Segment)
                segmentIndexMap[ seedSegmentIndex ] = length(segmentList1)
                # do not need to travase the subtree2 
                # so don't put children to seed list
            end 
        end 
    end 
    connectivityMatrix1 = spzeros(Bool, length(segmentList1), length(segmentList1))
    for (segmentIndexInSubtree1, segmentIndexInOriginalTree) in 
                                                    enumerate(segmentIndexListInSubtree1)
        parentSegmentIndexInOriginalTree = 
                                get_parent_segment_index(self, segmentIndexInOriginalTree)
        if parentSegmentIndexInOriginalTree > 0
            parentSegmentIndexInSubtree1 = segmentIndexMap[ parentSegmentIndexInOriginalTree ]
            connectivityMatrix1[ parentSegmentIndexInSubtree1, segmentIndexInSubtree1 ] = true 
        end 
    end 
    subtree1 = Neuron( segmentList1, connectivityMatrix1 )
    return subtree1, subtree2
end


function remove_segments!(self::Neuron, removeSegmentIndexList::Union{Vector,Set})
    self = remove_segments( self, removeSegmentIndexList )
end 

function remove_segments(self::Neuron, removeSegmentIndexList::Vector{Int})
    remove_segments(self, Set{Int}( removeSegmentIndexList ))
end 

function remove_segments(self::Neuron, removeSegmentIndexList::Set{Int})
    @assert !(1 in removeSegmentIndexList) "should not contain the root segment!"

    segmentList = get_segment_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    newSegmentList = Vector{Segment}()
    newSegmentIndexList = Vector{Int}()
    segmentIndexMap = Dict{Int,Int}()
    # the first one is the root segment
    seedSegmentIndexList = [1]
    # build segmentList1 without the splitout segment, will deal with it later
    while !isempty( seedSegmentIndexList )
        seedSegmentIndex = pop!( seedSegmentIndexList )
        if !(seedSegmentIndex in removeSegmentIndexList) 
            # should contain this segment
            push!(newSegmentList, segmentList[seedSegmentIndex])
            push!(newSegmentIndexList, seedSegmentIndex)
            # map the old segment index to the new segment index
            segmentIndexMap[ seedSegmentIndex ] = length(newSegmentList)
            childrenSegmentIndexList = get_children_segment_index_list(self, seedSegmentIndex)
            append!(seedSegmentIndexList, childrenSegmentIndexList)
        end 
    end 
    newConnectivityMatrix = spzeros(Bool, length(newSegmentList), length(newSegmentList))
    for (newSegmentIndex, segmentIndexInOriginalTree) in 
                                                    enumerate(newSegmentIndexList)
        parentSegmentIndexInOriginalTree = 
                                get_parent_segment_index(self, segmentIndexInOriginalTree)
        if parentSegmentIndexInOriginalTree > 0
            newParentSegmentIndex = segmentIndexMap[ parentSegmentIndexInOriginalTree ]
            newConnectivityMatrix[ newParentSegmentIndex, newSegmentIndex ] = true 
        end 
    end 
    return Neuron( newSegmentList, newConnectivityMatrix )
end 

"""
    remove_subtree_in_soma(self::Neuron)

remove the subtree which is inside the soma, which is an artifact of TEASAR algorithm
"""
function remove_subtree_in_soma( self::Neuron )
    removeSegmentIndexList = Vector{Int}()
    rootSegmentIndex = get_root_segment_index( self )
    rootNode = get_root_node( self )
    childrenSegmentIndexList = get_children_segment_index_list( self, rootSegmentIndex )
    for segmentIndex in childrenSegmentIndexList 
        terminalNodeList = get_terminal_node_list( self; startSegmentIndex=segmentIndex )
        if all( map(n->Segments.get_nodes_distance(n, rootNode) < rootNode[4]*2, 
                                                                terminalNodeList ) )
            println("remove segment: $segmentIndex")
            push!(removeSegmentIndexList, segmentIndex)
        end 
    end
    return remove_segments( self, removeSegmentIndexList )
end

function remove_hair( self::Neuron; radiiScale::Float32 = Float32(2),
                                    minDistance::Float32 = Float32(2000) )
    segmentList = get_segment_list(self)
    removeSegmentIndexList = Vector{Int}()
    for terminalSegmentIndex in get_terminal_segment_index_list(self)
        parentSegmentIndex = get_parent_segment_index(self, terminalSegmentIndex)
        terminalNode = segmentList[ terminalSegmentIndex ][end]
        parentNode = segmentList[ parentSegmentIndex ][end]
        distance = Segments.get_nodes_distance(terminalNode, parentNode)
        if distance < radiiScale*parentNode[4] || distance < minDistance
            println("remove segment $(terminalSegmentIndex) with distance of $(distance) and radius of $(parentNode[4])")
            push!(removeSegmentIndexList, terminalSegmentIndex)
        end 
    end 
    return remove_segments(self, removeSegmentIndexList)
end

"""
    remove_terminal_blobs( self::Neuron )
some terminal segmentation was fragmented. a lot of blobs was attatched to dendrite.
The blob terminal segment path length is normally smaller than the distance to parent dendrite.
"""
function remove_terminal_blobs( self::Neuron )
    terminalSegmentIndexList = get_terminal_segment_index_list( self )
    blobTerminalSegmentIndexList = Vector{Int}()
    for index in terminalSegmentIndexList 
        segment = self[index]
        distance2Parent = 0.0
        try 
            parentSegment = self[ get_parent_segment_index(self, index) ]
            distance2Parent = Segments.get_nodes_distance( segment[1], parentSegment[end] )
        catch err 
            @assert get_parent_segment_index(self, index) < 1
            warn("this terminal segment is root segment!")
        end 
        segmentInnerPathLength = Segments.get_path_length( segment )
        if distance2Parent > segmentInnerPathLength
            println("remove terminal blob. distance to parent: $(distance2Parent). segment inner path length: $(segmentInnerPathLength).")
            push!(blobTerminalSegmentIndexList, index)
        end 
    end 
    return remove_segments( self, blobTerminalSegmentIndexList )
end

"""
    remove_redundent_nodes( self::Neuron )
if neighboring node is the same, remove one of them 
"""
function remove_redundent_nodes(self::Neuron)
    newSegmentList = Vector{Segment}()
    sizehint!(newSegmentList, get_num_segments(self))
    removeSegmentIndexList = Vector{Int}()
    segmentList = get_segment_list(self)
    for (segmentIndex, segment) in enumerate( segmentList )
        Segments.remove_redundent_nodes!( segment )
        parentSegmentIndex = get_parent_segment_index(self, segmentIndex)
        if parentSegmentIndex > 0 
            parentSegment = segmentList[ parentSegmentIndex ]
            if segment[1] == parentSegment[end]
                segment = Segments.remove_node(segment, 1)
            end
        end 
        push!(newSegmentList, segment) 
        if length(segment) == 0
            # empty segment, should be removed later
            push!(removeSegmentIndexList, segmentIndex)
        end 
    end 
    newNeuron = Neuron( newSegmentList, get_connectivity_matrix(self) )
    return remove_segments(newNeuron, removeSegmentIndexList)
end 
########################## type convertion ####################
"""
    NodeNets.NodeNet( self::Neuron )
transform to NodeNet, the first node is the root node.
"""
function NodeNets.NodeNet(self::Neuron)
    nodeList = get_node_list( self )
    edges = get_edge_list( self )
    error("not fully implemented")
end 

function SWCs.SWC(self::Neuron)
    # initialize swc, which is a list of point objects
    swc = SWCs.SWC()
    # the node index of each segment, will be used for connecting the child segment
    segmentEndNodeIndexList = get_segment_end_node_index_list(self)
    # connectivity matrix of segments
    segmentConnectivityMatrix = get_connectivity_matrix(self)

    for (segmentIndex, segment) in enumerate( get_segment_list(self) )
        for (nodeIndex, node) in enumerate( Segments.get_node_list( segment ))
            # type, x, y, z, r, parent
            # in default, this was assumed that the connection is inside a segment, 
            # the parent is simply the previous node 
            pointObj = SWCs.PointObj( Segments.get_class(segment), node[1], node[2], node[3], node[4], length(swc) )
            
            if nodeIndex == 1
                # the first node should connect to other segment or be root node
                if length(swc) == 0
                    pointObj.parent = -1
                else
                    # find the connected parent segment index
                    parentSegmentIndex = get_parent_segment_index(self, segmentIndex)
                    # the node index of parent segment end
                    pointObj.parent= segmentEndNodeIndexList[ parentSegmentIndex ]
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
    # total number of bytes
    num_bytes = 4 + 4 + 4*3*get_node_num(self) + 4*2*length(get_edge_num(self))
    buffer = IOBuffer( num_bytes )
    # write the number of vertex
    write(buffer, UInt32(get_node_num(self)))
    write(buffer, UInt32(get_edge_num(self)))
    # write the node coordinates
    for node in get_node_list(self)
        write(buffer, [node[1:3]...])
    end 
    for edge in get_edges( self )
        write(buffer, UInt32( edge[1] ))
        write(buffer, UInt32( edge[2] ))
    end
    #bin = Vector{UInt8}(takebuf_array(buffer))
    bin = Vector{UInt8}(take!(buffer))
    close(buffer)
    return bin 
end 

############################### manipulations ##########################################
"""
    attach_pre_synapses!(self::Neuron, synapseTable::DataFrame)
make sure that the synapse table contains *only* the presynapses of this neuron 
"""
function attach_pre_synapses!(self::Neuron, synapseTable::DataFrame)
    for row in DataFrames.eachrow( synapseTable )
        synapse = Synapse(row)
        attach_pre_synapse!(self, synapse)
    end 
end 

"""
    attach_post_synapses!(self::Neuron, synapseTable::DataFrame)
make sure that the synapse table contains *only* the postsynapses of this neuron 
"""
function attach_post_synapses!(self::Neuron, synapseTable::DataFrame)
    for row in DataFrames.eachrow( synapseTable )
        synapse = Synapse(row)
        attach_post_synapse!(self, synapse)
    end 
end 


function attach_pre_synapse!(self::Neuron, synapse::Segments.Synapse)
    preSynapticCoordinate = Synapses.get_pre_synaptic_coordinate( synapse )
    segmentIndex, nodeIndexInSegment = find_closest_node_index( self, 
                                                                    preSynapticCoordinate )
    segmentList = get_segment_list(self)
    Segments.attach_pre_synapse!(segmentList[segmentIndex], nodeIndexInSegment, synapse) 
end 

function attach_post_synapse!(self::Neuron, synapse::Segments.Synapse)
    postSynapticCoordinate = Synapses.get_post_synaptic_coordinate( synapse )
    segmentIndex, nodeIndexInSegment = find_closest_node_index( self, 
                                                                postSynapticCoordinate )
    segmentList = get_segment_list(self)
    Segments.attach_post_synapse!(segmentList[segmentIndex], nodeIndexInSegment, synapse) 
end 


function downsample_nodes(self::Neuron; nodeNumStep::Int=DOWNSAMPLE_NODE_NUM_STEP) 
	@assert nodeNumStep > 1
    segmentList = get_segment_list(self)
    newSegmentList = Vector{Segment}()
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
    Neuron(newSegmentList, get_connectivity_matrix(self))
end

"""
    resample!( self::Neuron, resampleDistance::Float32 )
resampling the segments by a fixed point distance 
"""
function resample(self::Neuron, resampleDistance::Float32)
    newSegmentList = Vector{Segment}()
    segmentList = get_segment_list(self)
    local nodeList::Vector 
    for (index, segment) in enumerate(segmentList)
        parentSegmentIndex = get_parent_segment_index(self, index)
        if parentSegmentIndex > 0 
            parentNode = segmentList[parentSegmentIndex][end]
            nodeList = [parentNode, Segments.get_node_list(segment) ...]
        else 
            # no parent segment
            nodeList = Segments.get_node_list(segment) 
        end 
        nodeList = resample(nodeList, resampleDistance) 
        newSegment = Segment(nodeList[2:end]; class=Segments.get_class(segment))
        push!(newSegmentList, newSegment)
    end 
    Neuron(newSegmentList, get_connectivity_matrix(self))
end

function resample(nodeList::Vector{NTuple{4,Float32}}, resampleDistance::Float32)
    # always keep the first node
    ret = [nodeList[1]]

    # walk through the nodes
    walkedDistance = ZERO_FLOAT32 
    requiredDistance = resampleDistance  
    @inbounds for index in 1:length(nodeList)-1
        n1 = [nodeList[index  ] ...]
        n2 = [nodeList[index+1] ...]
        nodePairDistance = norm(n1[1:3].-n2[1:3])
        remainingDistance = nodePairDistance
        finishedDistance = ZERO_FLOAT32
        while remainingDistance > requiredDistance 
            # walk a step 
            # interperate the distance 
            interpretRatio = finishedDistance / nodePairDistance 
            interpretedNode = n1 .+ (n2.-n1) .* interpretRatio 
            push!(ret, (interpretedNode...))

            # update remaining distance 
            finishedDistance += requiredDistance   
            remainingDistance -= requiredDistance  
            requiredDistance = resampleDistance 
        end 
        requiredDistance = resampleDistance - remainingDistance 
    end 
    ret 
end

function find_closest_node_index(self::Neuron, coordinate::NTuple{3,Float32})
    segmentList = get_segment_list(self)
    
    # initialization 
    closestSegmentIndex = 0
    closestNodeIndexInSegment = 0
    distance = typemax(Float32)
    
    @assert !isempty(segmentList)
    
    for (segmentIndex, segment) in enumerate(segmentList)
        # distance from bounding box give the maximum bound of closest distance
        # this was used to filter out far away segments quickly
        # no need to compute all the distances between nodes
        # this is the lower bound of the distance
        bbox_distance = Segments.get_bounding_box_distance(segment, coordinate)
        if bbox_distance < distance 
            d, nodeIndexInSegment = Segments.distance_from(segment, coordinate )
            if d < distance 
                distance = d
                closestSegmentIndex = segmentIndex 
                closestNodeIndexInSegment = nodeIndexInSegment 
            end 
        end 
    end
    closestSegmentIndex, closestNodeIndexInSegment 
end 

function find_closest_terminal_node_index(self::Neuron, 
                                          nodeNet::NodeNet, 
                                          collectedFlagVec::Vector{Bool})
    # initialization 
    closestNodeIndex = 0
    closestTerminalSegmentIndex = 0
    closestDistance = typemax(Float32)

    terminalSegmentIndexList = get_terminal_segment_index_list(self)
    segmentList = get_segment_list(self)

    nodeList = NodeNets.get_node_list(nodeNet)
    connectivityMatrix = NodeNets.get_connectivity_matrix(nodeNet)
    @assert length(nodeList) == length(collectedFlagVec)

    for nodeIndex in 1:length(nodeList)
        # check all the uncollected terminal nodes 
        if !collectedFlagVec[nodeIndex] && NodeNets.is_terminal_node(nodeNet, nodeIndex)
            node = nodeList[nodeIndex]
            for terminalSegmentIndex in terminalSegmentIndexList
                terminalSegment = segmentList[terminalSegmentIndex]
                terminalNode = terminalSegment[end]
                # euclidean distance 
                distance = norm([map(-, terminalNode[1:3], node[1:3])...])
                if distance < closestDistance
                    closestDistance = distance
                    closestNodeIndex = nodeIndex
                    closestTerminalSegmentIndex = terminalSegmentIndex
                end
            end
        end
    end
    @assert closestNodeIndex > 0
    @assert closestDistance < typemax(Float32)
    return closestNodeIndex, closestTerminalSegmentIndex, closestDistance  
end 

function find_closest_node_index(neuron::Neuron, nodeNet::NodeNet, 
                                 collectedFlagVec::Vector{Bool})
    segmentList = get_segment_list(neuron)
    nodes = NodeNets.get_node_list(nodeNet)
    find_closest_node_index(segmentList, nodes, collectedFlagVec)
end 

"""
    find_closest_node_index(segmentList::Vector{Segment}, nodes::Vector{NTuple{4,Float32}})

find the uncollected node which is closest to the segment list 
"""
function find_closest_node_index(segmentList::Vector{Segment}, 
                                 nodes::Vector{NTuple{4,Float32}},
                                 collectedFlagVec ::Vector{Bool})
    # initialization 
    closestNodeIndex = 0
    closestSegmentIndex = 0
    closestNodeIndexInSegment = 0
    distance = typemax(Float32)
    
    @assert !isempty(segmentList)
    @assert length(nodes) == length(collectedFlagVec)
    @assert !all(collectedFlagVec)
    
    for (nodeIndex, node) in enumerate(nodes)
        if !collectedFlagVec[nodeIndex]
            for (segmentIndex, segment) in enumerate(segmentList)
                # distance from bounding box give the maximum bound of closest distance
                # this was used to filter out far away segments quickly
                # no need to compute all the distances between nodes
                # this is the lower bound of the distance
                bbox_distance = Segments.get_bounding_box_distance(segment, node)
                if bbox_distance < distance 
                    d, nodeIndexInSegment = Segments.distance_from(segment, node)
                    if d < distance 
                        distance = d
                        closestNodeIndex = nodeIndex 
                        closestSegmentIndex = segmentIndex 
                        closestNodeIndexInSegment = nodeIndexInSegment 
                    end 
                end 
            end
        end 
    end
    # the segmentIndex should be inside the segmentList
    #@assert closestSegmentIndex > 0
    #@assert closestSegmentIndex <= length(segmentList)
    return closestNodeIndex, closestSegmentIndex, closestNodeIndexInSegment 
end 

function add_offset(self::Neuron, offset::Union{Vector, Tuple})
    segmentList = Vector{Segment}()
    for segment in self.segmentList 
        newSegment = Segments.add_offset(segment, offset)
        push!(segmentList, newSegment)
    end 
    Neuron(segmentList, self.connectivityMatrix)
end 

end # module
