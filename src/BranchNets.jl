module BranchNets
include("Branches.jl")
using .Branches 
using ..TEASAR.NodeNets
using ..TEASAR.SWCs

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
    numBranchPoint = NodeNets.get_num_branch_point(nodeNet)
    # number of branches, assumes that this is acyclic graph
    numBranches = numBranchPoint + 2
    # initializae the net
    branchList = Vector{Branch}()
    connectivityMatrix = spzeros(Bool, numBranches, numBranches)
    net = BranchNet(branchList, connectivityMatrix)

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
    while !all(collectedFlagVec)
        # there exist some uncollected nodes 
        # find the uncollected node that is nearest to the branch list as seed 
        nearestNodeIndex, nearestBranchIndex, nearestNodeIndexInBranch = 
                            find_nearest_node_index(branchList, nodes, collectedFlagVec)
        if isempty(seedNodeIndexList)
            # this is not the first time, the main net is already built
            seedNodeIndexList = [nearestNodeIndex]
        end  
        while !isempty(seedNodeIndexList)
            # grow the net from seed until there is branching point
            seedNodeIndex = pop!(seedNodeIndexList)
            # grow a subnet
            subnet = BranchNet!(seedNodeIndex, nodeNet, collectedFlagVec)
            # merge the subnet to main net
            net = merge!(net, subnet, nearestBranchIndex, nearestNodeIndexInBranch)
        end 
    end 
    net     
end 

function BranchNet!(seedNodeIndex::Integer, nodeNet::NodeNet, collectedFlagVec::Vector{Bool})
    # initialization
    branchList = Vector{Branch}()
    parentBranchIndexList = Vector{Int}()
    childBranchIndexList  = Vector{Int}()

    nodes = NodeNets.get_node_list(nodeNet)
    nodesConnectivityMatrix = NodeNets.get_connectivity_matrix( nodeNet )
    
    # the seed of branch should record both seed node index and the the parent branch index
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
            @show countnz(collectedFlagVec)
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
    connectivityMatrix = sparse(parentBranchIndexList, childBranchIndexList, 
                                ones(Bool,length(childBranchIndexList)))
    BranchNet(branchList, connectivityMatrix)
end 

"""
    get_num_branches(self::BranchNet)

get the number of branches  
"""
function get_num_branches(self::BranchNet) length(self.branchList) end 
function get_branch_list(self::BranchNet) self.branchList end 
function get_connectivity_matrix(self::BranchNet) self.connectivityMatrix end 
"""
    get_branch_length_list( self::BranchNet )

get a vector of Integer, which represent the length of each branch 
"""
function get_branch_length_list( self::BranchNet ) map(length, get_branch_list(self)) end

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
merge two nets at a specific branch location
"""
function Base.merge!(self::BranchNet, other::BranchNet, nearestBranchIndex::Integer, 
                     nearestNodeIndexInBranch::Integer)
    if isempty(self) 
        return other  
    elseif nearestNodeIndexInBranch == length(self.branchList[nearestBranchIndex])
        # the connection point is a branching point, no need to break branch
        # just stitch the new branches and rebuild the connection matrix
        mergedBranchList = vcat( self.branchList, other.branchList )
        # total number of branches
        num_branches1 = get_num_branches(self)
        num_branches2 = get_num_branches(other)
        total_num_branches = get_num_branches(self) + get_num_branches(other)
        mergedConnectivityMatrix = spzeros(Bool, total_num_branches, total_num_branches)
        mergedConnectivityMatrix[1:num_branches1, 1:num_branches1] = 
                                                        self.connectivityMatrix 
        mergedConnectivityMatrix[num_branches1+1:end, num_branches1+1:end] = 
                                                        other.connectivityMatrix
        # connect the two nets, the row index is the parent and the column is the child
        mergedConnectivityMatrix[nearestBranchIndex, num_branches1+1] = true
        return BranchNet(mergedBranchList, mergedConnectivityMatrix)
    else 
        # need to break the nearest branch and rebuild connectivity matrix
        num_branches1 = get_num_branches(self)
        num_branches2 = get_num_branches(other)
        total_num_branches = num_branches1 + num_branches2 + 1
        mergedConnectivityMatrix = spzeros(Bool, total_num_branches, total_num_branches)
        
        # need to break the branch and then stitch the new branches
        branch1, branch2 = split(self.branchList[nearestBranchIndex], 
                                                nearestNodeIndexInBranch)
        self.branchList[nearestBranchIndex] = branch1 
        mergedConnectivityMatrix[1:num_branches1, 1:num_branches1] = 
                                                                self.connectivityMatrix
        # reconnect the breaked two branches
        push!(self.branchList, branch2)
        mergedConnectivityMatrix[nearestBranchIndex, num_branches1+1] = true 

        # merge the other net
        push!(self.branchList, other.branchList)
        mergedConnectivityMatrix[num_branches1+2:end, num_branches1+2:end] = 
                                                                other.connectivityMatrix 
        
        # create new merged net
        return BranchNet(self.branchList, mergedConnectivityMatrix)
    end
end 

function Base.isempty(self::BranchNet)    isempty(self.branchList) end 

########################## type convertion ####################
function SWCs.SWC(self::BranchNet)
    # initialize swc, which is a list of point objects
    swc = SWC()
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
                    parentBranchIndexList,_ = findnz(branchConnectivityMatrix[:, branchIndex])
                    # always have and only have one parent branch
                    @assert length(parentBranchIndexList) == 1
                    parentBranchIndex = parentBranchIndexList[1]
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
    find_nearest_node_index(branchList::Vector{Branch}, nodes::Vector{NTuple{4,Float32}})

find the node which is nearest to the branch list 
"""
function find_nearest_node_index(branchList::Vector{Branch}, 
                                 nodes::Vector{NTuple{4,Float32}},
                                 collectedFlagVec ::Vector{Bool})
    ret = (0,0,0)
    distance = typemax(Float32)
    @assert length(nodes) == length(collectedFlagVec)
    for (nodeIndex, node) in enumerate(nodes)
        if !collectedFlagVec[nodeIndex]
            for (branchIndex, branch) in enumerate(branchList)
                # distance from bounding box give the maximum bound of nearest distance
                # this was used to filter out far away branches quickly
                # no need to compute all the distances between nodes
                bbox_distance = Branches.get_bounding_box_distance(branch, node)
                if bbox_distance < distance 
                    d, nodeIndexInBranch = Branches.distance_from(branch, node)
                    if d < distance 
                        distance = d
                        ret = (nodeIndex, branchIndex, nodeIndexInBranch)
                    end 
                end 
            end
        end 
    end 
    ret 
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
