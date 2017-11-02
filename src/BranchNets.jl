module BranchNets
include("Branches.jl")
using .Branches 
using ..RealNeuralNetworks.NodeNets
using ..RealNeuralNetworks.SWCs

const ONE_UINT32 = UInt32(1)
const EXPANSION = (ONE_UINT32, ONE_UINT32, ONE_UINT32)

export BranchNet

type BranchNet 
    # x,y,z,r
    branchList ::Vector{Branch}
    connectivityMatrix ::SparseMatrixCSC{Bool, Int}
end 

"""
    BranchNet
a neuron modeled by interconnected branches 
"""
function BranchNet(nodeNet::NodeNet)
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
    branchNet = BranchNet!(rootNodeIndex, nodeNet, collectedFlagVec)
    #return branchNet # only use the main branch for evaluation

    while !all(collectedFlagVec)
        # there exist some uncollected nodes 
        # find the uncollected node that is nearest to the branch list as seed
        branchList = get_branch_list( branchNet )
        
        # there still exist uncollected nodes
        seedNodeIndex, nearestBranchIndex, nearestNodeIndexInBranch = 
                    find_nearest_node_index(branchList, nodes, collectedFlagVec)
        @assert nearestBranchIndex <= length(branchList) 
        # grow a new net from seed, and mark the collected nodes 
        subnet = BranchNet!(seedNodeIndex, nodeNet, collectedFlagVec) 
        # merge the subnet to main net 
        branchNet = merge( branchNet, subnet, 
                                nearestBranchIndex, nearestNodeIndexInBranch)
    end 
    branchNet     
end 

"""
    BranchNet!(seedNodeIndex::Integer, nodeNet::NodeNet, collectedFlagVec::Vector{Bool})

build a net from a seed node using connected component
mark the nodes in this new net as collected, so the collectedFlagVec was changed.
"""
function BranchNet!(seedNodeIndex::Integer, nodeNet::NodeNet, 
                                            collectedFlagVec::Vector{Bool})
    # initialization
    branchList = Vector{Branch}()
    parentBranchIndexList = Vector{Int}()
    childBranchIndexList  = Vector{Int}()

    nodes = NodeNets.get_node_list(nodeNet)
    nodesConnectivityMatrix = NodeNets.get_connectivity_matrix( nodeNet )
    
    # the seed of branch should record both seed node index 
    # and the the parent branch index
    # the parent branch index of root branch is -1 
    branchSeedList = [(seedNodeIndex, -1)]
    # depth first search
    while !isempty(branchSeedList)
        seedNodeIndex, branchParentIndex = pop!(branchSeedList)
        # grow a branch from seed
        nodeListInBranch = Vector{NTuple{4, Float32}}()
        seedNodeIndexList = [seedNodeIndex]
        while true
            # construct this branch
            seedNodeIndex = pop!(seedNodeIndexList)
            # push the seed node
            push!(nodeListInBranch, nodes[seedNodeIndex])
            # label this node index as collected
            collectedFlagVec[ seedNodeIndex ] = true 

            # find the connected nodes
            connectedNodeIndexList,_ = findnz(nodesConnectivityMatrix[:, seedNodeIndex])
            # exclude the collected nodes
            connectedNodeIndexList = connectedNodeIndexList[ !collectedFlagVec[connectedNodeIndexList] ] 

            if length(connectedNodeIndexList) == 1
                # belong to the same branch
                push!(nodeListInBranch, nodes[ connectedNodeIndexList[1] ])
                push!(seedNodeIndexList, connectedNodeIndexList[1])
            else
                # terminal branching point or multiple branching points
                # finish constructing this branch
                branch = Branch(nodeListInBranch)
                push!(branchList, branch)
                if branchParentIndex != -1
                    # this is not the root branch, establish the branch connection
                    push!(parentBranchIndexList, branchParentIndex)
                    push!(childBranchIndexList,  length(branchList))
                end 
                # seed new branches
                # if this is terminal branch, no seed will be pushed
                for index in connectedNodeIndexList 
                    push!(branchSeedList, (index, length(branchList)))
                end 
                break
            end 
        end 
    end
    @assert length(parentBranchIndexList) == length(childBranchIndexList)
    # note that the connectivity matrix should be a square matrix for easier use
    connectivityMatrix = spzeros(Bool, length(branchList), length(branchList))
    tempConnectivityMatrix = sparse(parentBranchIndexList, childBranchIndexList, 
                                ones(Bool,length(childBranchIndexList)))
    connectivityMatrix[1:size(tempConnectivityMatrix,1), 
                       1:size(tempConnectivityMatrix,2)] = tempConnectivityMatrix 
    BranchNet(branchList, connectivityMatrix)
end 

function BranchNet{T}( seg::Array{T,3}; obj_id::T = convert(T,1), 
                        expansion::NTuple{3,UInt32}=EXPANSION )
    nodeNet = NodeNet( seg; obj_id = obj_id, expansion = expansion )
    BranchNet( nodeNet )
end

function BranchNet( swc::SWC )
    nodeNet = NodeNet( swc )
    BranchNet( nodeNet )
end 

"""
    BranchNet( swcString::AbstractString )
transform string from swc content to BranchNet 
"""
function BranchNet( swcString::AbstractString )
    swc = SWC(swcString)
    BranchNet( swc )
end

######################### IO ################
function load_swc( fileName::AbstractString )
    swcString = readstring( fileName )
    BranchNet( swcString )
end

function load_gzip_swc( fileName::AbstractString )
    swc = SWCs.load_gzip_swc( fileName )
    BranchNet( swc )
end 

function save(self::BranchNet, fileName::AbstractString)
    swc = SWC(self)
    SWCs.save( swc, fileName )
end

function save_gzip_swc( self::BranchNet, fileName::AbstractString )
    SWCs.save_gzip_swc( SWCs.SWC(self), fileName )
end 

####################### properties ##############

"""
    get_root_branch_index( self::BranchNet )

the first branch should be the root branch, and the first node should be the root node 
"""
function get_root_branch_index(self::BranchNet) 1 end 
function get_root_branch(self::BranchNet) 
    get_branch_list(self)[get_root_branch_index(self)]
end 

"""
    get_root_node( self::BranchNet )
the returned root node is a tuple of Float64, which represents x,y,z,r
"""
function get_root_node( self::BranchNet ) 
    rootBranchIndex = get_root_branch_index(self)
    rootBranch = get_branch_list(self)[ rootBranchIndex ]
    rootBranch[1]
end 

"""
    get_num_branches(self::BranchNet)
"""
function get_num_branches(self::BranchNet) length(self.branchList) end

"""
    get_num_branching_points(self::BranchNet)
"""
function get_num_branching_points(self::BranchNet)
    numBranchingPoint = 0
    for index in 1:get_num_branches(self)
        childrenBranchIndexList = get_children_branch_index_list(self, index)
        if length(childrenBranchIndexList) > 0
            numBranchingPoint += 1
        end 
    end 
    numBranchingPoint
end 

function get_branch_list(self::BranchNet) self.branchList end 
function get_connectivity_matrix(self::BranchNet) self.connectivityMatrix end 

"""
    get_branch_order_list( self::BranchNet )
following the Julia indexing style, the root branch order is 1.
"""
function get_branch_order_list(self::BranchNet)
    branchOrderList = Vector{Int}()
    sizehint!(branchOrderList, get_num_branches(self))
    
    index2order = Dict{Int,Int}()
    # the first one is branch index, the second one is branch order 
    seedBranchIndexOrderList = Vector{NTuple{2,Int}}([(1,1)])
    while !isempty( seedBranchIndexOrderList )
        branchIndex, branchOrder = pop!(seedBranchIndexOrderList)
        index2order[branchIndex] = branchOrder 
        childrenBranchIndexList = get_children_branch_index_list(self, branchIndex)
        for childBranchIndex in childrenBranchIndexList 
            push!(seedBranchIndexOrderList, (childBranchIndex, branchOrder+1))
        end 
    end 
    
    for index in 1:get_num_branches( self )
        push!(branchOrderList, index2order[index])
    end 
    branchOrderList
end 

"""
    get_branch_length_list( self::BranchNet )

get a vector of Integer, which represent the length of each branch 
"""
function get_branch_length_list( self::BranchNet ) map(length, get_branch_list(self)) end

"""
    get_node_num( self::BranchNet )
get total number of nodes  
"""
function get_node_num( self::BranchNet )
    sum(map(length, get_branch_list(self)))
end 

"""
    get_node_list(self::BranchNet)
get the node list. the first one is root node.
"""
function get_node_list( self::BranchNet )
    nodeList = Vector{NTuple{4, Float32}}()
    sizehint!(nodeList, get_node_num(self))
    for branch in get_branch_list(self)
        append!(nodeList, Branches.get_node_list(branch))
    end 
    nodeList
end 

"""
    get_edge_list( self::BranchNet )
get the edges with type of Vector{NTuple{2, Int}}
"""
function get_edge_list( self::BranchNet )
    edgeList = Vector{NTuple{2,Int}}()
    branchStartNodeIndexList = Vector{Int}()
    branchStopNodeIndexList = Vector{Int}()
    # total number of nodes
    nodeNum = 0
    for (branchIndex, branch) in enumerate(get_branch_list(self))
        push!(branchStartNodeIndexList, nodeNum+1)
        # build the edges inside branch
        for nodeIndex in nodeNum+1:nodeNum+length(branch)-1
            push!(edgeList, (nodeIndex, nodeIndex+1))
        end 
        # update node number
        nodeNum += length(branch)
        push!(branchStopNodeIndexList,  nodeNum)
    end 
    # add branch connections
    parentBranchIndexList, childBranchIndexList, _ = findnz( get_connectivity_matrix(self) )
    for (index, parentBranchIndex) in enumerate( parentBranchIndexList )
        childBranchIndex = childBranchIndexList[ index ]
        parentNodeIndex = branchStopNodeIndexList[ parentBranchIndex ]
        childNodeIndex  = branchStartNodeIndexList[ childBranchIndex ]
        push!( edgeList, (parentNodeIndex, childNodeIndex) )
    end 
    edgeList 
end 

function get_path_to_root_length(self::BranchNet, branchIndex::Integer)
    path2RootLength = 0.0
    branchPathLengthList = get_branch_path_length_list( self )
    while true 
        path2RootLength += branchPathLengthList[ branchIndex ]
        parentBranchIndex = get_parent_branch_index(self, branchIndex )
        if parentBranchIndex < 1 
            # root branch do not have parent 
            break 
        else
            branchIndex = parentBranchIndex 
        end 
    end
    path2RootLength 
end

"""
    get_branch_path_length_list(self::BranchNet)
get euclidean path length of each branch 
"""
function get_branch_path_length_list( self::BranchNet )
    ret = Vector{Float64}()
    for (index, branch) in enumerate( get_branch_list(self) )
        branchPathLength = Branches.get_path_length( branch )
        # add the edge length to parent node
        parentBranchIndex = get_parent_branch_index(self, index)
        if parentBranchIndex > 0
            parentBranch = get_branch_list(self)[ parentBranchIndex ]
            parentNode = parentBranch[ end ]
            node = branch[1]
            branchPathLength += norm( [map((x,y)->x-y, node[1:3], parentNode[1:3])...] )
        end 
        push!(ret, branchPathLength)
    end 
    ret 
end 

"""
    get_branch_end_node_index_list( self::BranchNet )

get a vector of integer, which represent the node index of branch end 
"""
function get_branch_end_node_index_list( self::BranchNet ) cumsum( get_branch_length_list(self) ) end 

"""
    get_num_nodes(self::BranchNet)

get number of nodes 
"""
function get_num_nodes( self::BranchNet )
    branchList = get_branch_list(self)
    sum(map(length, branchList))
end 

"""
    get_mass_center( self::BranchNet )
mass center was computed simply as center of all nodes 
"""
function get_mass_center( self::BranchNet )
    nodeList = get_node_list(self)
    x = sum(map(n->n[1], nodeList)) / length(nodeList)
    y = sum(map(n->n[2], nodeList)) / length(nodeList)
    z = sum(map(n->n[3], nodeList)) / length(nodeList)
    (x,y,z)
end

"""
    get_asymmetry( self::BranchNet )
asymmetry was measured by the euclidean distance between root node and mass center 
"""
function get_asymmetry( self::BranchNet )
    massCenter = get_mass_center( self )
    root = get_root_node( self )
    norm( [massCenter...] .- [root[1:3]...] )
end 

"""
    get_children_branch_index_list(self::BranchNet, parentBranchIndex::Integer)
"""
function get_children_branch_index_list(self::BranchNet, parentBranchIndex::Integer)
    childrenBranchIndexList,_ = findnz(self.connectivityMatrix[parentBranchIndex, :])
    childrenBranchIndexList 
end

function get_parent_branch_index( self::BranchNet, childBranchIndex::Integer )
    parentBranchIndexList,_ = findnz(self.connectivityMatrix[:, childBranchIndex])
    @assert length(parentBranchIndexList) <= 1
    if isempty( parentBranchIndexList ) 
        # no parent, this is a root branch
        return 0
    else 
        return parentBranchIndexList[1]
    end 
end

"""
    get_subtree_branch_index_list( self, branchInde )
get the branch index list of subtree 
"""
function get_subtree_branch_index_list( self::BranchNet, branchIndex::Integer )
    @assert branchIndex > 0 && branchIndex <= get_num_branches(self)
    subtreeBranchIndexList = Vector{Integer}()
    seedBranchIndexList = [branchIndex]
    while !isempty( seedBranchIndexList )
        seedBranchIndex = pop!( seedBranchIndexList )
        push!(subtreeBranchIndexList, seedBranchIndex)
        childrenBranchIndexList = get_children_branch_index_list( self, seedBranchIndex )
        append!(seedBranchIndexList, childrenBranchIndexList)
    end 
    subtreeBranchIndexList 
end 

function get_terminal_branch_index_list( self::BranchNet; startBranchIndex::Integer = 1 )
    terminalBranchIndexList = Vector{Int}()
    branchIndexList = [ startBranchIndex ]
    while !isempty( branchIndexList )
        branchIndex = pop!( branchIndexList )
        childrenBranchIndexList = get_children_branch_index_list( self, branchIndex )
        if isempty( childrenBranchIndexList ) 
            # this is a terminal branch
            push!(terminalBranchIndexList, branchIndex )
        else
            append!(branchIndexList, childrenBranchIndexList)
        end 
    end 
    terminalBranchIndexList 
end

function get_terminal_node_list( self::BranchNet; startBranchIndex::Integer = 1 )
    terminalBranchIndexList = get_terminal_branch_index_list( self )
    map( x -> Branches.get_node_list(x)[end], get_branch_list(self)[ terminalBranchIndexList ] )
end 

"""
    get_branching_angle( self::BranchNet, branchIndex::Integer; nodeDistance::AbstractFloat  )
if the node is too close the angle might not be accurate. For example, if they are voxel neighbors, the angle will alwasys be 90 or 45 degree. Thus, we have nodeDistance here to help make the nodes farther away.
Note that the acos returns angle in the format of radiens.
"""
function get_branching_angle( self::BranchNet, branchIndex::Integer; nodeDistance::Real = 5000.0 )
    branch = self[branchIndex]
    parentBranch = self[ get_parent_branch_index(self, branchIndex) ]
    branchingNode = parentBranch[end]
    parentNode = parentBranch[end-1]
    for node in parentBranch 
        if Branches.get_nodes_distance(node, branchingNode) < nodeDistance
            parentNode = node 
            break 
        end 
    end
    childNode = branch[1]
    for index in length(branch):1
        node = branch[index]
        if Branches.get_nodes_distance(node, branchingNode) < nodeDistance 
            childNode = node 
            break 
        end 
    end 
    # compute the angle among the three nodes using the definition of dot product
    # get the two vectors
    v1 = map(-, parentNode[1:3], branchingNode[1:3] )
    v2 = map(-, branchingNode[1:3], childNode[1:3] )
    # normalize the vector
    nv1 = normalize([v1...]) 
    nv2 = normalize([v2...])
    return acos( dot(nv1, nv2) )
end 

"""
    get_sholl_number(self::BranchNet, shollRadius::AbstractFloat)
suppose there is a sphere centered on the root node. The sholl number is the number of contact points of the neuron and sphere. 
"""
function get_sholl_number(self::BranchNet, shollRadius::AbstractFloat)
    shollNum = 0
    root = get_root_node(self)
    nodeList = get_node_list( self )
    edgeList = get_edge_list( self )
    for edge in edgeList 
        node1 = nodeList[ edge[1] ]
        node2 = nodeList[ edge[2] ]
        if  Branches.get_nodes_distance(root, node1) > shollRadius != 
            Branches.get_nodes_distance(root, node2) > shollRadius 
            shollNum += 1
        end 
    end 
    shollNum 
end 

function get_sholl_number_list(self::BranchNet, shollRadiusList::Vector)
    shollNumList = zeros(Float64, length(shollRadiusList))
    for (index, shollRadius) in enumerate(shollRadiusList)
        shollNumList = get_sholl_number(self, shollRadius) 
    end 
    shollNumList 
end

"""
    get_sholl_number_list(self::BranchNet, radiusStep::AbstractFloat; numStep::Integer=0)
if the number of steps is 0, we'll auto compute the step number according to the farthest terminal node.
"""
function get_sholl_number_list(self::BranchNet, radiusStep::Real;
                               numStep::Integer=0)
    shollNumList = Vector{Float64}()
    if numStep==0
        # automatically computing step number 
        root = get_root_node(self)
        terminalNodeList = get_terminal_node_list(self)
        terminal2RootDistanceList = map(n->Branches.get_nodes_distance(root,n), terminalNodeList)
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


############################### Base functions ###################################
function Base.getindex(self::BranchNet, index::Integer)
    get_branch_list(self)[index]
end 
"""
merge two nets at a specific branch location.
the root node of second net will be connected to the first net
"""
function Base.merge(self::BranchNet, other::BranchNet, 
                    nearestBranchIndex::Integer, nearestNodeIndexInBranch::Integer)
    @assert !isempty(self)
    @assert !isempty(other)
    branchList1 = get_branch_list( self  )
    branchList2 = get_branch_list( other )
    num_branches1 = get_num_branches(self)
    num_branches2 = get_num_branches(other)
    
    if nearestNodeIndexInBranch == length(branchList1[nearestBranchIndex])
        # the connection point is the end of a branch, no need to break branch
        # merge the root branch of second net to the first net
        mergedBranchList1 = copy(branchList1)
        # assume that the first branch is the root branch
        mergedBranchList1[nearestBranchIndex] = 
                                merge(mergedBranchList[nearestBranchIndex], branchList2[1])
        mergedBranchList = vcat( mergedBranchList1, branchList2[2:end] )

        # total number of branches
        total_num_branches = num_branches1 + num_branches2 - 1 
        
        mergedConnectivityMatrix = spzeros(Bool, total_num_branches, total_num_branches)
        mergedConnectivityMatrix[1:size(self.connectivityMatrix,1), 
                                 1:size(self.connectivityMatrix,2)] = 
                                                                self.connectivityMatrix 
        # do not include the connection of root in net2
        mergedConnectivityMatrix[
                num_branches1+1+1 : num_branches1+1+size(other.connectivityMatrix,1)-1, 
                num_branches1+1 : num_branches1+1+size(other.connectivityMatrix,2)-1] = 
                                                        other.connectivityMatrix[2:end, :]
        # reestablish the connection of root2
        childrenBranchIndexList2 = get_children_branch_index_list(other, 1)
        for childBranchIndex2 in childrenBranchIndexList2
            mergedConnectivityMatrix[ nearestBranchIndex, 
                                      num_branches1+childBranchIndex-1 ] = true 
        end 
        return BranchNet(mergedBranchList, mergedConnectivityMatrix)
    else 
        # need to break the nearest branch and rebuild connectivity matrix
        total_num_branches = num_branches1 + 1 + num_branches2 
        mergedBranchList = branchList1 
        mergedConnectivityMatrix = spzeros(Bool, total_num_branches, total_num_branches)
        
        # need to break the branch and then stitch the new branches
        branchPart1, branchPart2 = split(branchList1[nearestBranchIndex], 
                                                    nearestNodeIndexInBranch)
        mergedBranchList[nearestBranchIndex] = branchPart1 
        mergedConnectivityMatrix[1:size(self.connectivityMatrix,1), 
                                 1:size(self.connectivityMatrix,2)] = 
                                                        self.connectivityMatrix 
        # reconnect the breaked two branches
        push!(mergedBranchList, branchPart2)
        mergedConnectivityMatrix[nearestBranchIndex, num_branches1+1] = true 

        # redirect the children branches to branchPart2
        childrenBranchIndexList = get_children_branch_index_list(self, nearestBranchIndex)
        for childBranchIndex in childrenBranchIndexList
            # remove old connection
            mergedConnectivityMatrix[nearestBranchIndex, childBranchIndex] = false
            # build new connection
            mergedConnectivityMatrix[num_branches1+1, childBranchIndex] = true 
        end 

        # merge the other net
        mergedBranchList = vcat(mergedBranchList, branchList2)
        mergedConnectivityMatrix[
                num_branches1+2 : num_branches1+2+size(other.connectivityMatrix,1)-1, 
                num_branches1+2 : num_branches1+2+size(other.connectivityMatrix,2)-1] = 
                                                                other.connectivityMatrix

        # establish the connection between two nets
        mergedConnectivityMatrix[nearestBranchIndex, num_branches1+2] = true

        # create new merged net
        return BranchNet(mergedBranchList, mergedConnectivityMatrix)
    end
end 

function Base.isempty(self::BranchNet)    isempty(self.branchList) end 

"""
    Base.split(self::BranchNet, splitBranchIndex::Integer; nodeIndexInBranch::Integer=0)

split the branch net into two branchNets
the nodeIndexInBranch will be included in the first main net including the original root node  
"""
function Base.split(self::BranchNet, splitBranchIndex::Integer; nodeIndexInBranch::Integer=1)
    if splitBranchIndex == 1 && nodeIndexInBranch < 1 
        return BranchNet(), self
    end 

    branchList = get_branch_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    branchIndexListInSubtree2 = get_subtree_branch_index_list( self, splitBranchIndex )
    subtree1Branch, subtree2RootBranch = split(branchList[splitBranchIndex], 
                                               nodeIndexInBranch)

    # build the split out subtree2, which do not contain the root node
    branchList2 = [subtree2RootBranch]
    append!(branchList2, branchList[ branchIndexListInSubtree2 ])
    branchIndexMap = Dict{Int,Int}()
    sizehint!(branchIndexMap, length(branchIndexListInSubtree2))
    for (branchIndexInSubtree2, mainTreeBranchIndex ) in 
                                                    enumerate(branchIndexListInSubtree2)
        branchIndexMap[ mainTreeBranchIndex ] = branchIndexInSubtree2 
    end
    connectivityMatrix2 = spzeros( Bool, length(branchList2), length(branchList2) )
    for (branchIndexInSubtree2, branchIndexInOriginalTree) in 
                                            enumerate(branchIndexListInSubtree2)
        parentBranchIndexInOriginalTree = 
                                get_parent_branch_index(self, branchIndexInOriginalTree)
        if haskey(branchIndexMap, parentBranchIndexInOriginalTree)
            parentBranchIndexInSubtree2 = branchIndexMap[ parentBranchIndexInOriginalTree ]
            connectivityMatrix2[ parentBranchIndexInSubtree2, branchIndexInSubtree2 ] = true 
        else 
            # this is the root branch of new tree 
            println("find root of new tree in the original tree : $(parentBranchIndexInOriginalTree)")
        end 
    end 
    subtree2 = BranchNet( branchList2, connectivityMatrix2 )

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    branchList1 = Vector{Branch}()
    branchIndexListInSubtree1 = Vector{Int}()
    branchIndexMap = Dict{Int,Int}()
    # the first one is the root branch
    seedBranchIndexList = [1]
    # build branchList1 without the splitout branch, will deal with it later
    while !isempty( seedBranchIndexList )
        seedBranchIndex = pop!( seedBranchIndexList )
        if seedBranchIndex != splitBranchIndex 
            push!(branchList1, branchList[seedBranchIndex])
            push!(branchIndexListInSubtree1, seedBranchIndex)
            # map the old branch index to the new branch index
            # 3 => 2 means the 3rd branch of old tree is the 2nd branch of new tree
            branchIndexMap[ seedBranchIndex ] = length(branchList1)
            childrenBranchIndexList = get_children_branch_index_list(self, seedBranchIndex)
            append!(seedBranchIndexList, childrenBranchIndexList)
        else 
            # deal with the split out branch
            if nodeIndexInBranch != 1
                @assert length(subtree1Branch) != 0
                push!(branchList1, subtree1Branch)
                branchIndexMap[ seedBranchIndex ] = length(branchList1)
                # do not need to travase the subtree2 
                # so don't put children to seed list
            end 
        end 
    end 
    connectivityMatrix1 = spzeros(Bool, length(branchList1), length(branchList1))
    for (branchIndexInSubtree1, branchIndexInOriginalTree) in 
                                                    enumerate(branchIndexListInSubtree1)
        parentBranchIndexInOriginalTree = 
                                get_parent_branch_index(self, branchIndexInOriginalTree)
        if parentBranchIndexInOriginalTree > 0
            parentBranchIndexInSubtree1 = branchIndexMap[ parentBranchIndexInOriginalTree ]
            connectivityMatrix1[ parentBranchIndexInSubtree1, branchIndexInSubtree1 ] = true 
        end 
    end 
    subtree1 = BranchNet( branchList1, connectivityMatrix1 )
    return subtree1, subtree2
end


function remove_branches!(self::BranchNet, removeBranchIndexList::Union{Vector,Set})
    self = remove_branches( self, removeBranchIndexList )
end 

function remove_branches(self::BranchNet, removeBranchIndexList::Vector{Int})
    remove_branches(self, Set{Int}( removeBranchIndexList ))
end 

function remove_branches(self::BranchNet, removeBranchIndexList::Set{Int})
    @assert !(1 in removeBranchIndexList) "should not contain the root branch!"

    branchList = get_branch_list(self)
    connectivityMatrix = get_connectivity_matrix(self)

    # rebuild the main tree without the subtree2. 
    # the main tree still contains the old root
    newBranchList = Vector{Branch}()
    newBranchIndexList = Vector{Int}()
    branchIndexMap = Dict{Int,Int}()
    # the first one is the root branch
    seedBranchIndexList = [1]
    # build branchList1 without the splitout branch, will deal with it later
    while !isempty( seedBranchIndexList )
        seedBranchIndex = pop!( seedBranchIndexList )
        if !(seedBranchIndex in removeBranchIndexList) 
            # should contain this branch
            push!(newBranchList, branchList[seedBranchIndex])
            push!(newBranchIndexList, seedBranchIndex)
            # map the old branch index to the new branch index
            branchIndexMap[ seedBranchIndex ] = length(newBranchList)
            childrenBranchIndexList = get_children_branch_index_list(self, seedBranchIndex)
            append!(seedBranchIndexList, childrenBranchIndexList)
        end 
    end 
    newConnectivityMatrix = spzeros(Bool, length(newBranchList), length(newBranchList))
    for (newBranchIndex, branchIndexInOriginalTree) in 
                                                    enumerate(newBranchIndexList)
        parentBranchIndexInOriginalTree = 
                                get_parent_branch_index(self, branchIndexInOriginalTree)
        if parentBranchIndexInOriginalTree > 0
            newParentBranchIndex = branchIndexMap[ parentBranchIndexInOriginalTree ]
            newConnectivityMatrix[ newParentBranchIndex, newBranchIndex ] = true 
        end 
    end 
    return BranchNet( newBranchList, newConnectivityMatrix )
end 

"""
    remove_subtree_in_soma(self::BranchNet)

remove the subtree which is inside the soma, which is an artifact of TEASAR algorithm
"""
function remove_subtree_in_soma( self::BranchNet )
    removeBranchIndexList = Vector{Int}()
    rootBranchIndex = get_root_branch_index( self )
    rootNode = get_root_node( self )
    childrenBranchIndexList = get_children_branch_index_list( self, rootBranchIndex )
    for branchIndex in childrenBranchIndexList 
        terminalNodeList = get_terminal_node_list( self; startBranchIndex=branchIndex )
        if all( map(n->Branches.get_nodes_distance(n, rootNode) < rootNode[4]*2, 
                                                                terminalNodeList ) )
            println("remove branch: $branchIndex")
            push!(removeBranchIndexList, branchIndex)
        end 
    end
    return remove_branches( self, removeBranchIndexList )
end

function remove_hair( self::BranchNet )
    branchList = get_branch_list(self)
    removeBranchIndexList = Vector{Int}()
    for terminalBranchIndex in get_terminal_branch_index_list(self)
        parentBranchIndex = get_parent_branch_index(self, terminalBranchIndex)
        terminalNode = branchList[ terminalBranchIndex ][end]
        parentNode = branchList[ parentBranchIndex ][end]
        distance = Branches.get_nodes_distance(terminalNode, parentNode)
        if distance < 2*parentNode[4]
            println("remove branch $(terminalBranchIndex) with distance of $(distance) and radius of $(parentNode[4])")
            push!(removeBranchIndexList, terminalBranchIndex)
        end 
    end 
    return remove_branches(self, removeBranchIndexList)
end

########################## type convertion ####################
"""
    NodeNets.NodeNet( self::BranchNet )
transform to NodeNet, the first node is the root node.
"""
function NodeNets.NodeNet(self::BranchNet)
    nodeList = get_node_list( self )
    edges = get_edge_list( self )
    error("not fully implemented")
end 

function SWCs.SWC(self::BranchNet)
    # initialize swc, which is a list of point objects
    swc = SWCs.SWC()
    # the node index of each branch, will be used for connecting the child branch
    branchEndNodeIndexList = get_branch_end_node_index_list(self)
    # connectivity matrix of branches
    branchConnectivityMatrix = get_connectivity_matrix(self)

    for (branchIndex, branch) in enumerate( get_branch_list(self) )
        for (nodeIndex, node) in enumerate( Branches.get_node_list( branch ))
            # type, x, y, z, r, parent
            # in default, this was assumed that the connection is inside a branch, 
            # the parent is simply the previous node 
            pointObj = SWCs.PointObj( Branches.get_class(branch), node[1], node[2], node[3], node[4], length(swc) )
            
            if nodeIndex == 1
                # the first node should connect to other branch or be root node
                if length(swc) == 0
                    pointObj.parent = -1
                else
                    # find the connected parent branch index
                    parentBranchIndex = get_parent_branch_index(self, branchIndex)
                    # the node index of parent branch end
                    pointObj.parent= branchEndNodeIndexList[ parentBranchIndex ]
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
function get_neuroglancer_precomputed(self::BranchNet)
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

############################### TEASAR algorithm functions #############################
"""
    find_nearest_node_index(branchList::Vector{Branch}, nodes::Vector{NTuple{4,Float32}})

find the uncollected node which is nearest to the branch list 
"""
function find_nearest_node_index(branchList::Vector{Branch}, 
                                 nodes::Vector{NTuple{4,Float32}},
                                 collectedFlagVec ::Vector{Bool})
    # initialization 
    nearestNodeIndex = 0
    nearestBranchIndex = 0
    nearestNodeIndexInBranch = 0
    distance = typemax(Float32)
    
    @assert !isempty(branchList)
    @assert length(nodes) == length(collectedFlagVec)
    @assert !all(collectedFlagVec)
    
    for (nodeIndex, node) in enumerate(nodes)
        if !collectedFlagVec[nodeIndex]
            for (branchIndex, branch) in enumerate(branchList)
                # distance from bounding box give the maximum bound of nearest distance
                # this was used to filter out far away branches quickly
                # no need to compute all the distances between nodes
                # this is the lower bound of the distance
                bbox_distance = Branches.get_bounding_box_distance(branch, node)
                if bbox_distance < distance 
                    d, nodeIndexInBranch = Branches.distance_from(branch, node)
                    if d < distance 
                        distance = d
                        nearestNodeIndex = nodeIndex 
                        nearestBranchIndex = branchIndex 
                        nearestNodeIndexInBranch = nodeIndexInBranch 
                    end 
                end 
            end
        end 
    end
    # the branchIndex should be inside the branchList
    @assert nearestBranchIndex > 0
    @assert nearestBranchIndex <= length(branchList)
    return nearestNodeIndex, nearestBranchIndex, nearestNodeIndexInBranch 
end 

function add_offset(self::BranchNet, offset::Union{Vector, Tuple})
    branchList = Vector{Branch}()
    for branch in self.branchList 
        newBranch = Branches.add_offset(branch, offset)
        push!(branchList, newBranch)
    end 
    BranchNet(branchList, self.connectivityMatrix)
end 

end # module
