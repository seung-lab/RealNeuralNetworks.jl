module Neurons
include("Segments.jl")
using .Segments
using ..RealNeuralNetworks.NodeNets
using ..RealNeuralNetworks.SWCs
using LsqFit
using ImageFiltering
using OffsetArrays

const ONE_UINT32 = UInt32(1)
const ZERO_FLOAT32 = Float32(0)
const EXPANSION = (ONE_UINT32, ONE_UINT32, ONE_UINT32)

export Neuron

mutable struct Neuron 
    # x,y,z,r, the coordinates should be in physical coordinates
    segmentList ::Vector{Branch}
    connectivityMatrix ::SparseMatrixCSC{Bool, Int}
end 

"""
    Neuron
a neuron modeled by interconnected segmentes 
"""
function Neuron(nodeNet::NodeNet)
#    @save "nodenet.jld2" nodeNet 
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
        # find the uncollected node that is nearest to the segment list as seed
        segmentList = get_segment_list( neuron )
        
        # there still exist uncollected nodes
        seedNodeIndex, nearestBranchIndex, nearestNodeIndexInBranch = 
                    find_nearest_node_index(segmentList, nodes, collectedFlagVec)
        @assert nearestBranchIndex <= length(segmentList) 
        # grow a new net from seed, and mark the collected nodes 
        subnet = Neuron!(seedNodeIndex, nodeNet, collectedFlagVec) 
        # merge the subnet to main net 
        neuron = merge( neuron, subnet, 
                                nearestBranchIndex, nearestNodeIndexInBranch)
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
    segmentList = Vector{Branch}()
    parentBranchIndexList = Vector{Int}()
    childBranchIndexList  = Vector{Int}()

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
        nodeListInBranch = Vector{NTuple{4, Float32}}()
        seedNodeIndexList = [seedNodeIndex]
        while true
            # construct this segment
            seedNodeIndex = pop!(seedNodeIndexList)
            # push the seed node
            push!(nodeListInBranch, nodes[seedNodeIndex])
            # label this node index as collected
            collectedFlagVec[ seedNodeIndex ] = true 

            # find the connected nodes
            connectedNodeIndexList,_ = findnz(nodesConnectivityMatrix[:, seedNodeIndex])
            # exclude the collected nodes
            connectedNodeIndexList = connectedNodeIndexList[ .!collectedFlagVec[connectedNodeIndexList] ] 

            if length(connectedNodeIndexList) == 1
                # belong to the same segment
                push!(nodeListInBranch, nodes[ connectedNodeIndexList[1] ])
                push!(seedNodeIndexList, connectedNodeIndexList[1])
            else
                # terminal segmenting point or multiple segmenting points
                # finish constructing this segment
                segment = Branch(nodeListInBranch)
                push!(segmentList, segment)
                if segmentParentIndex != -1
                    # this is not the root segment, establish the segment connection
                    push!(parentBranchIndexList, segmentParentIndex)
                    push!(childBranchIndexList,  length(segmentList))
                end 
                # seed new segmentes
                # if this is terminal segment, no seed will be pushed
                for index in connectedNodeIndexList 
                    push!(segmentSeedList, (index, length(segmentList)))
                end 
                break
            end 
        end 
    end
    @assert length(parentBranchIndexList) == length(childBranchIndexList)
    # note that the connectivity matrix should be a square matrix for easier use
    connectivityMatrix = spzeros(Bool, length(segmentList), length(segmentList))
    tempConnectivityMatrix = sparse(parentBranchIndexList, childBranchIndexList, 
                                ones(Bool,length(childBranchIndexList)))
    connectivityMatrix[1:size(tempConnectivityMatrix,1), 
                       1:size(tempConnectivityMatrix,2)] = tempConnectivityMatrix 
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
    rootBranchIndex = get_root_segment_index(self)
    rootBranch = get_segment_list(self)[ rootBranchIndex ]
    rootBranch[1]
end 

"""
    get_num_segmentes(self::Neuron)
"""
function get_num_segmentes(self::Neuron) length(self.segmentList) end

"""
    get_num_segmenting_points(self::Neuron)
"""
function get_num_segmenting_points(self::Neuron)
    numBranchingPoint = 0
    for index in 1:get_num_segmentes(self)
        childrenBranchIndexList = get_children_segment_index_list(self, index)
        if length(childrenBranchIndexList) > 0
            numBranchingPoint += 1
        end 
    end 
    numBranchingPoint
end 

function get_segment_list(self::Neuron) self.segmentList end 
function get_connectivity_matrix(self::Neuron) self.connectivityMatrix end 

"""
    get_segment_order_list( self::Neuron )
following the Julia indexing style, the root segment order is 1.
"""
function get_segment_order_list(self::Neuron)
    segmentOrderList = Vector{Int}()
    sizehint!(segmentOrderList, get_num_segmentes(self))
    
    index2order = Dict{Int,Int}()
    # the first one is segment index, the second one is segment order 
    seedBranchIndexOrderList = Vector{NTuple{2,Int}}([(1,1)])
    while !isempty( seedBranchIndexOrderList )
        segmentIndex, segmentOrder = pop!(seedBranchIndexOrderList)
        index2order[segmentIndex] = segmentOrder 
        childrenBranchIndexList = get_children_segment_index_list(self, segmentIndex)
        for childBranchIndex in childrenBranchIndexList 
            push!(seedBranchIndexOrderList, (childBranchIndex, segmentOrder+1))
        end 
    end 
    
    for index in 1:get_num_segmentes( self )
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
    parentBranchIndexList, childBranchIndexList, _ = findnz( get_connectivity_matrix(self) )
    for (index, parentBranchIndex) in enumerate( parentBranchIndexList )
        childBranchIndex = childBranchIndexList[ index ]
        parentNodeIndex = segmentStopNodeIndexList[ parentBranchIndex ]
        childNodeIndex  = segmentStartNodeIndexList[ childBranchIndex ]
        push!( edgeList, (parentNodeIndex, childNodeIndex) )
    end 
    edgeList 
end 

function get_path_to_root_length(self::Neuron, segmentIndex::Integer; segmentPathLengthList::Vector = [])
    path2RootLength = 0.0
    if isempty(segmentPathLengthList)
        segmentPathLengthList = get_segment_path_length_list( self )
    end 
    while true 
        path2RootLength += segmentPathLengthList[ segmentIndex ]
        parentBranchIndex = get_parent_segment_index(self, segmentIndex )
        if parentBranchIndex < 1 
            # root segment do not have parent 
            break 
        else
            segmentIndex = parentBranchIndex 
        end 
    end
    path2RootLength 
end

"""
    get_segment_path_length_list(self::Neuron)
get euclidean path length of each segment 
"""
function get_segment_path_length_list( self::Neuron )
    ret = Vector{Float64}()
    for (index, segment) in enumerate( get_segment_list(self) )
        segmentPathLength = Segments.get_path_length( segment )
        # add the edge length to parent node
        parentBranchIndex = get_parent_segment_index(self, index)
        if parentBranchIndex > 0
            parentBranch = get_segment_list(self)[ parentBranchIndex ]
            parentNode = parentBranch[ end ]
            node = segment[1]
            segmentPathLength += norm( [map((x,y)->x-y, node[1:3], parentNode[1:3])...] )
        end 
        push!(ret, segmentPathLength)
    end 
    ret 
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
    get_children_segment_index_list(self::Neuron, parentBranchIndex::Integer)
"""
function get_children_segment_index_list(self::Neuron, parentBranchIndex::Integer)
    childrenBranchIndexList,_ = findnz(self.connectivityMatrix[parentBranchIndex, :])
    childrenBranchIndexList 
end

function get_parent_segment_index( self::Neuron, childBranchIndex::Integer )
    parentBranchIndexList,_ = findnz(self.connectivityMatrix[:, childBranchIndex])
    @assert length(parentBranchIndexList) <= 1
    if isempty( parentBranchIndexList ) 
        # no parent, this is a root segment
        return 0
    else 
        return parentBranchIndexList[1]
    end 
end

"""
    get_subtree_segment_index_list( self, segmentInde )
get the segment index list of subtree 
"""
function get_subtree_segment_index_list( self::Neuron, segmentIndex::Integer )
    @assert segmentIndex > 0 && segmentIndex <= get_num_segmentes(self)
    subtreeBranchIndexList = Vector{Integer}()
    seedBranchIndexList = [segmentIndex]
    while !isempty( seedBranchIndexList )
        seedBranchIndex = pop!( seedBranchIndexList )
        push!(subtreeBranchIndexList, seedBranchIndex)
        childrenBranchIndexList = get_children_segment_index_list( self, seedBranchIndex )
        append!(seedBranchIndexList, childrenBranchIndexList)
    end 
    subtreeBranchIndexList 
end 

function get_terminal_segment_index_list( self::Neuron; startBranchIndex::Integer = 1 )
    terminalBranchIndexList = Vector{Int}()
    segmentIndexList = [ startBranchIndex ]
    while !isempty( segmentIndexList )
        segmentIndex = pop!( segmentIndexList )
        childrenBranchIndexList = get_children_segment_index_list( self, segmentIndex )
        if isempty( childrenBranchIndexList ) 
            # this is a terminal segment
            push!(terminalBranchIndexList, segmentIndex )
        else
            append!(segmentIndexList, childrenBranchIndexList)
        end 
    end 
    terminalBranchIndexList 
end

function get_terminal_node_list( self::Neuron; startBranchIndex::Integer = 1 )
    terminalBranchIndexList = get_terminal_segment_index_list( self )
    map( x -> Segments.get_node_list(x)[end], get_segment_list(self)[ terminalBranchIndexList ] )
end 

"""
    get_segmenting_angle( self::Neuron, segmentIndex::Integer; nodeDistance::AbstractFloat  )
if the node is too close the angle might not be accurate. For example, if they are voxel neighbors, the angle will alwasys be 90 or 45 degree. Thus, we have nodeDistance here to help make the nodes farther away.
Note that the acos returns angle in the format of radiens.
"""
function get_segmenting_angle( self::Neuron, segmentIndex::Integer; nodeDistance::Real = 5000.0 )
    segment = self[segmentIndex]
    parentBranchIndex = get_parent_segment_index(self, segmentIndex)
    if parentBranchIndex < 1
        # this is root segment
        return 0.0
    end
    parentBranch = self[ get_parent_segment_index(self, segmentIndex) ]
    if length(parentBranch) == 1 || length(segment) == 1
        return 0.0
    end 
    segmentingNode = parentBranch[end]
    parentNode = parentBranch[end-1]
    for node in parentBranch 
        if Segments.get_nodes_distance(node, segmentingNode) < nodeDistance
            parentNode = node 
            break 
        end 
    end
    if parentNode == segmentingNode
        warn("parent node is the same with segmenting node: $(segmentingNode)")
        return 0.0
    end 

    childNode = segment[1]
    for index in length(segment):1
        node = segment[index]
        if Segments.get_nodes_distance(node, segmentingNode) < nodeDistance 
            childNode = node 
            break 
        end 
    end 
    if childNode == segmentingNode 
        warn("child node is the same with segmenting node: $(segmentingNode)")
        return 0.0 
    end 

    # compute the angle among the three nodes using the definition of dot product
    # get the two vectors
    v1 = map(-, parentNode[1:3], segmentingNode[1:3] )
    v2 = map(-, segmentingNode[1:3], childNode[1:3] )
    # normalize the vector
    nv1 = normalize([v1...]) 
    nv2 = normalize([v2...])
    #@show nv1, nv2
    #@show childNode, segmentingNode, parentNode
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
"""
function get_mask(self::Neuron, voxelSize::Union{Tuple, Vector})
    nodeList = get_node_list(self)
    voxelList = Vector{NTuple{3, Int}}()
    for node in nodeList 
        voxelCoordinate = map((x,y)->round(Int, x/y), node[1:3], voxelSize)
        push!(voxelList, (voxelCoordinate...))
    end 
    boundingBox = Segments.BoundingBox( voxelList )
    range = Segments.BoundingBoxes.get_unit_range(boundingBox)
    sz = size(boundingBox)
    # initialize the map 
    mask = zeros(Float64, sz)
    mask = OffsetArray(mask, range...)

    @unsafe for voxel in voxelList
        # add a small number to make sure that there is no 0
        mask[voxel...] = true 
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
    @time densityMap = imfilter(densityMap, kernel)
    # normalize by total path length of the neuron
    totalPathLength = get_total_path_length(self)
    const scale = totalPathLength / norm(densityMap[:]) 
    densityMap .*= scale 
    densityMap
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
                    nearestBranchIndex::Integer, nearestNodeIndexInBranch::Integer)
    @assert !isempty(self)
    @assert !isempty(other)
    segmentList1 = get_segment_list( self  )
    segmentList2 = get_segment_list( other )
    num_segmentes1 = get_num_segmentes(self)
    num_segmentes2 = get_num_segmentes(other)
    
    if nearestNodeIndexInBranch == length(segmentList1[nearestBranchIndex])
        println("connecting to a segment end...")
        childrenBranchIndexList = get_children_segment_index_list(self, nearestBranchIndex)
        if length(childrenBranchIndexList) > 0
            println("the nearest segment have children, do not merge the root segment")
            total_num_segmentes = num_segmentes1 + num_segmentes2 
            mergedBranchList = vcat(segmentList1, segmentList2)
            @assert length(mergedBranchList) == total_num_segmentes 
            mergedConnectivityMatrix = 
                        spzeros(Bool, total_num_segmentes, total_num_segmentes)
            mergedConnectivityMatrix[   1:size(self.connectivityMatrix,1), 
                                        1:size(self.connectivityMatrix,2)] = 
                                                                self.connectivityMatrix 
            # do not include the connection of root in net2
            mergedConnectivityMatrix[num_segmentes1+1 : end, 
                                     num_segmentes1+1 : end] = other.connectivityMatrix 
            mergedConnectivityMatrix[nearestBranchIndex, num_segmentes1+1] = true
            return Neuron(mergedBranchList, mergedConnectivityMatrix)
        else 
            println("the nearest segment is a terminal segment, merge the root segment of 'other'")
            if num_segmentes2 == 1
                #  only one segment of the other net
                mergedBranchList = copy(segmentList1)
                mergedBranchList[nearestBranchIndex] = 
                            merge(mergedBranchList[nearestBranchIndex], segmentList2[1])
                return Neuron(mergedBranchList, self.connectivityMatrix)
            else 
                # the connection point is the end of a segment, no need to break segment
                # merge the root segment of second net to the first net
                # assume that the first segment is the root segment
                mergedBranchList = vcat( segmentList1, segmentList2[2:end] )
                mergedBranchList[nearestBranchIndex] = 
                            merge(mergedBranchList[nearestBranchIndex], segmentList2[1])

                # total number of segmentes
                total_num_segmentes = num_segmentes1 + num_segmentes2 - 1 
                @assert length(mergedBranchList) == total_num_segmentes 
                
                mergedConnectivityMatrix = 
                                spzeros(Bool, total_num_segmentes, total_num_segmentes)
                mergedConnectivityMatrix[
                                    1:size(self.connectivityMatrix,1), 
                                    1:size(self.connectivityMatrix,2)] = 
                                                        self.connectivityMatrix 
                # do not include the connection of root in net2
                mergedConnectivityMatrix[num_segmentes1+1 : end, 
                                         num_segmentes1+1 : end] = 
                                                other.connectivityMatrix[2:end, 2:end]
                # reestablish the connection of root2
                childrenBranchIndexList2 = get_children_segment_index_list(other, 1)
                for childBranchIndex2 in childrenBranchIndexList2
                    mergedConnectivityMatrix[ nearestBranchIndex, 
                                              num_segmentes1+childBranchIndex2-1 ] = true 
                end
                return Neuron(mergedBranchList, mergedConnectivityMatrix)
            end 
        end 
    else 
        #println("need to break the nearest segment and rebuild connectivity matrix")
        total_num_segmentes = num_segmentes1 + 1 + num_segmentes2 
        mergedBranchList = segmentList1 
        mergedConnectivityMatrix = spzeros(Bool, total_num_segmentes, total_num_segmentes)
        
        # need to break the segment and then stitch the new segmentes
        segmentPart1, segmentPart2 = split(segmentList1[nearestBranchIndex], 
                                                    nearestNodeIndexInBranch)
        mergedBranchList[nearestBranchIndex] = segmentPart1 
        mergedConnectivityMatrix[1:size(self.connectivityMatrix,1), 
                                 1:size(self.connectivityMatrix,2)] = 
                                                        self.connectivityMatrix 
        # reconnect the breaked two segmentes
        push!(mergedBranchList, segmentPart2)
        mergedConnectivityMatrix[nearestBranchIndex, num_segmentes1+1] = true 

        # redirect the children segmentes to segmentPart2
        childrenBranchIndexList = get_children_segment_index_list(self, nearestBranchIndex)
        for childBranchIndex in childrenBranchIndexList
            # remove old connection
            mergedConnectivityMatrix[nearestBranchIndex, childBranchIndex] = false
            # build new connection
            mergedConnectivityMatrix[num_segmentes1+1, childBranchIndex] = true 
        end 

        # merge the other net
        mergedBranchList = vcat(mergedBranchList, segmentList2)
        mergedConnectivityMatrix[
                num_segmentes1+2 : num_segmentes1+2+size(other.connectivityMatrix,1)-1, 
                num_segmentes1+2 : num_segmentes1+2+size(other.connectivityMatrix,2)-1] = 
                                                                other.connectivityMatrix

        # establish the connection between two nets
        mergedConnectivityMatrix[nearestBranchIndex, num_segmentes1+2] = true

        # create new merged net
        return Neuron(mergedBranchList, mergedConnectivityMatrix)
    end
end 

function Base.isempty(self::Neuron)    isempty(self.segmentList) end 

"""
    Base.split(self::Neuron, splitBranchIndex::Integer; nodeIndexInBranch::Integer=0)

split the segment net into two neurons
the nodeIndexInBranch will be included in the first main net including the original root node  
"""
function Base.split(self::Neuron, splitBranchIndex::Integer; nodeIndexInBranch::Integer=1)
    if splitBranchIndex == 1 && nodeIndexInBranch < 1 
        return Neuron(), self
    end 

    segmentList = get_segment_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    segmentIndexListInSubtree2 = get_subtree_segment_index_list( self, splitBranchIndex )
    subtree1Branch, subtree2RootBranch = split(segmentList[splitBranchIndex], 
                                                                nodeIndexInBranch)

    # build the split out subtree2, which do not contain the root node
    segmentList2 = [subtree2RootBranch]
    append!(segmentList2, segmentList[ segmentIndexListInSubtree2 ])
    segmentIndexMap = Dict{Int,Int}()
    sizehint!(segmentIndexMap, length(segmentIndexListInSubtree2))
    for (segmentIndexInSubtree2, mainTreeBranchIndex ) in 
                                                    enumerate(segmentIndexListInSubtree2)
        segmentIndexMap[ mainTreeBranchIndex ] = segmentIndexInSubtree2 
    end
    connectivityMatrix2 = spzeros( Bool, length(segmentList2), length(segmentList2) )
    for (segmentIndexInSubtree2, segmentIndexInOriginalTree) in 
                                            enumerate(segmentIndexListInSubtree2)
        parentBranchIndexInOriginalTree = 
                                get_parent_segment_index(self, segmentIndexInOriginalTree)
        if haskey(segmentIndexMap, parentBranchIndexInOriginalTree)
            parentBranchIndexInSubtree2 = segmentIndexMap[ parentBranchIndexInOriginalTree ]
            connectivityMatrix2[ parentBranchIndexInSubtree2, segmentIndexInSubtree2 ] = true 
        else 
            # this is the root segment of new tree 
            println("find root of new tree in the original tree : $(parentBranchIndexInOriginalTree)")
        end 
    end 
    subtree2 = Neuron( segmentList2, connectivityMatrix2 )

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    segmentList1 = Vector{Branch}()
    segmentIndexListInSubtree1 = Vector{Int}()
    segmentIndexMap = Dict{Int,Int}()
    # the first one is the root segment
    seedBranchIndexList = [1]
    # build segmentList1 without the splitout segment, will deal with it later
    while !isempty( seedBranchIndexList )
        seedBranchIndex = pop!( seedBranchIndexList )
        if seedBranchIndex != splitBranchIndex 
            push!(segmentList1, segmentList[seedBranchIndex])
            push!(segmentIndexListInSubtree1, seedBranchIndex)
            # map the old segment index to the new segment index
            # 3 => 2 means the 3rd segment of old tree is the 2nd segment of new tree
            segmentIndexMap[ seedBranchIndex ] = length(segmentList1)
            childrenBranchIndexList = get_children_segment_index_list(self, seedBranchIndex)
            append!(seedBranchIndexList, childrenBranchIndexList)
        else 
            # deal with the split out segment
            if nodeIndexInBranch != 1
                @assert length(subtree1Branch) != 0
                push!(segmentList1, subtree1Branch)
                segmentIndexMap[ seedBranchIndex ] = length(segmentList1)
                # do not need to travase the subtree2 
                # so don't put children to seed list
            end 
        end 
    end 
    connectivityMatrix1 = spzeros(Bool, length(segmentList1), length(segmentList1))
    for (segmentIndexInSubtree1, segmentIndexInOriginalTree) in 
                                                    enumerate(segmentIndexListInSubtree1)
        parentBranchIndexInOriginalTree = 
                                get_parent_segment_index(self, segmentIndexInOriginalTree)
        if parentBranchIndexInOriginalTree > 0
            parentBranchIndexInSubtree1 = segmentIndexMap[ parentBranchIndexInOriginalTree ]
            connectivityMatrix1[ parentBranchIndexInSubtree1, segmentIndexInSubtree1 ] = true 
        end 
    end 
    subtree1 = Neuron( segmentList1, connectivityMatrix1 )
    return subtree1, subtree2
end


function remove_segmentes!(self::Neuron, removeBranchIndexList::Union{Vector,Set})
    self = remove_segmentes( self, removeBranchIndexList )
end 

function remove_segmentes(self::Neuron, removeBranchIndexList::Vector{Int})
    remove_segmentes(self, Set{Int}( removeBranchIndexList ))
end 

function remove_segmentes(self::Neuron, removeBranchIndexList::Set{Int})
    @assert !(1 in removeBranchIndexList) "should not contain the root segment!"

    segmentList = get_segment_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    newBranchList = Vector{Branch}()
    newBranchIndexList = Vector{Int}()
    segmentIndexMap = Dict{Int,Int}()
    # the first one is the root segment
    seedBranchIndexList = [1]
    # build segmentList1 without the splitout segment, will deal with it later
    while !isempty( seedBranchIndexList )
        seedBranchIndex = pop!( seedBranchIndexList )
        if !(seedBranchIndex in removeBranchIndexList) 
            # should contain this segment
            push!(newBranchList, segmentList[seedBranchIndex])
            push!(newBranchIndexList, seedBranchIndex)
            # map the old segment index to the new segment index
            segmentIndexMap[ seedBranchIndex ] = length(newBranchList)
            childrenBranchIndexList = get_children_segment_index_list(self, seedBranchIndex)
            append!(seedBranchIndexList, childrenBranchIndexList)
        end 
    end 
    newConnectivityMatrix = spzeros(Bool, length(newBranchList), length(newBranchList))
    for (newBranchIndex, segmentIndexInOriginalTree) in 
                                                    enumerate(newBranchIndexList)
        parentBranchIndexInOriginalTree = 
                                get_parent_segment_index(self, segmentIndexInOriginalTree)
        if parentBranchIndexInOriginalTree > 0
            newParentBranchIndex = segmentIndexMap[ parentBranchIndexInOriginalTree ]
            newConnectivityMatrix[ newParentBranchIndex, newBranchIndex ] = true 
        end 
    end 
    return Neuron( newBranchList, newConnectivityMatrix )
end 

"""
    remove_subtree_in_soma(self::Neuron)

remove the subtree which is inside the soma, which is an artifact of TEASAR algorithm
"""
function remove_subtree_in_soma( self::Neuron )
    removeBranchIndexList = Vector{Int}()
    rootBranchIndex = get_root_segment_index( self )
    rootNode = get_root_node( self )
    childrenBranchIndexList = get_children_segment_index_list( self, rootBranchIndex )
    for segmentIndex in childrenBranchIndexList 
        terminalNodeList = get_terminal_node_list( self; startBranchIndex=segmentIndex )
        if all( map(n->Segments.get_nodes_distance(n, rootNode) < rootNode[4]*2, 
                                                                terminalNodeList ) )
            println("remove segment: $segmentIndex")
            push!(removeBranchIndexList, segmentIndex)
        end 
    end
    return remove_segmentes( self, removeBranchIndexList )
end

function remove_hair( self::Neuron )
    segmentList = get_segment_list(self)
    removeBranchIndexList = Vector{Int}()
    for terminalBranchIndex in get_terminal_segment_index_list(self)
        parentBranchIndex = get_parent_segment_index(self, terminalBranchIndex)
        terminalNode = segmentList[ terminalBranchIndex ][end]
        parentNode = segmentList[ parentBranchIndex ][end]
        distance = Segments.get_nodes_distance(terminalNode, parentNode)
        if distance < 2*parentNode[4] || distance < 2000
            println("remove segment $(terminalBranchIndex) with distance of $(distance) and radius of $(parentNode[4])")
            push!(removeBranchIndexList, terminalBranchIndex)
        end 
    end 
    return remove_segmentes(self, removeBranchIndexList)
end

"""
    remove_terminal_blobs( self::Neuron )
some terminal segmentation was fragmented. a lot of blobs was attatched to dendrite.
The blob terminal segment path length is normally smaller than the distance to parent dendrite.
"""
function remove_terminal_blobs( self::Neuron )
    terminalBranchIndexList = get_terminal_segment_index_list( self )
    blobTerminalBranchIndexList = Vector{Int}()
    for index in terminalBranchIndexList 
        segment = self[index]
        distance2Parent = 0.0
        try 
            parentBranch = self[ get_parent_segment_index(self, index) ]
            distance2Parent = Segments.get_nodes_distance( segment[1], parentBranch[end] )
        catch err 
            @assert get_parent_segment_index(self, index) < 1
            warn("this terminal segment is root segment!")
        end 
        segmentInnerPathLength = Segments.get_path_length( segment )
        if distance2Parent > segmentInnerPathLength
            println("remove terminal blob. distance to parent: $(distance2Parent). segment inner path length: $(segmentInnerPathLength).")
            push!(blobTerminalBranchIndexList, index)
        end 
    end 
    return remove_segmentes( self, blobTerminalBranchIndexList )
end

"""
    remove_redundent_nodes( self::Neuron )
if neighboring node is the same, remove one of them 
"""
function remove_redundent_nodes(self::Neuron)
    newBranchList = Vector{Branch}()
    sizehint!(newBranchList, get_num_segmentes(self))
    removeBranchIndexList = Vector{Int}()
    segmentList = get_segment_list(self)
    for (segmentIndex, segment) in enumerate( segmentList )
        Segments.remove_redundent_nodes!( segment )
        parentBranchIndex = get_parent_segment_index(self, segmentIndex)
        if parentBranchIndex > 0 
            parentBranch = segmentList[ parentBranchIndex ]
            if segment[1] == parentBranch[end]
                segment = Segments.remove_node(segment, 1)
            end
        end 
        push!(newBranchList, segment) 
        if length(segment) == 0
            # empty segment, should be removed later
            push!(removeBranchIndexList, segmentIndex)
        end 
    end 
    newNeuron = Neuron( newBranchList, get_connectivity_matrix(self) )
    return remove_segmentes(newNeuron, removeBranchIndexList)
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
    # connectivity matrix of segmentes
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
                    parentBranchIndex = get_parent_segment_index(self, segmentIndex)
                    # the node index of parent segment end
                    pointObj.parent= segmentEndNodeIndexList[ parentBranchIndex ]
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
    resample!( self::Neuron, resampleDistance::Float32 )
resampling the segmentes by a fixed point distance 
"""
function resample(self::Neuron, resampleDistance::Float32)
    newBranchList = Vector{Branch}()
    segmentList = get_segment_list(self)
    local nodeList::Vector 
    for (index, segment) in enumerate(segmentList)
        parentBranchIndex = get_parent_segment_index(self, index)
        if parentBranchIndex > 0 
            parentNode = segmentList[parentBranchIndex][end]
            nodeList = [parentNode, Segments.get_node_list(segment) ...]
        else 
            # no parent segment
            nodeList = Segments.get_node_list(segment) 
        end 
        nodeList = resample(nodeList, resampleDistance) 
        newBranch = Branch(nodeList[2:end]; class=Segments.get_class(segment))
        push!(newBranchList, newBranch)
    end 
    Neuron(newBranchList, get_connectivity_matrix(self))
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

"""
    find_nearest_node_index(segmentList::Vector{Branch}, nodes::Vector{NTuple{4,Float32}})

find the uncollected node which is nearest to the segment list 
"""
function find_nearest_node_index(segmentList::Vector{Branch}, 
                                 nodes::Vector{NTuple{4,Float32}},
                                 collectedFlagVec ::Vector{Bool})
    # initialization 
    nearestNodeIndex = 0
    nearestBranchIndex = 0
    nearestNodeIndexInBranch = 0
    distance = typemax(Float32)
    
    @assert !isempty(segmentList)
    @assert length(nodes) == length(collectedFlagVec)
    @assert !all(collectedFlagVec)
    
    for (nodeIndex, node) in enumerate(nodes)
        if !collectedFlagVec[nodeIndex]
            for (segmentIndex, segment) in enumerate(segmentList)
                # distance from bounding box give the maximum bound of nearest distance
                # this was used to filter out far away segmentes quickly
                # no need to compute all the distances between nodes
                # this is the lower bound of the distance
                bbox_distance = Segments.get_bounding_box_distance(segment, node)
                if bbox_distance < distance 
                    d, nodeIndexInBranch = Segments.distance_from(segment, node)
                    if d < distance 
                        distance = d
                        nearestNodeIndex = nodeIndex 
                        nearestBranchIndex = segmentIndex 
                        nearestNodeIndexInBranch = nodeIndexInBranch 
                    end 
                end 
            end
        end 
    end
    # the segmentIndex should be inside the segmentList
    #@assert nearestBranchIndex > 0
    #@assert nearestBranchIndex <= length(segmentList)
    #@show distance
    #if distance > 1000 
    #    @show nearestNodeIndex, nearestBranchIndex, nearestNodeIndexInBranch
    #end 
    return nearestNodeIndex, nearestBranchIndex, nearestNodeIndexInBranch 
end 

function add_offset(self::Neuron, offset::Union{Vector, Tuple})
    segmentList = Vector{Branch}()
    for segment in self.segmentList 
        newBranch = Segments.add_offset(segment, offset)
        push!(segmentList, newBranch)
    end 
    Neuron(segmentList, self.connectivityMatrix)
end 

end # module
