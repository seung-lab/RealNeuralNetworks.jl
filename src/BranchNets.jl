module BranchNets
include("Branches.jl")
using .Branches 
using ..TEASAR.NodeNets

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
    nodes = NodeNets.get_nodes(nodeNet)
    radii = NodeNets.get_radii(nodeNet)
    # flags labeling whether this node was collected to the net
    collectedFlagVec = zeros(Bool, length(nodes))
    # connectivity matrix of nodes in nodeNet 
    nodesConnectivityMatrix  = NodeNets.get_connectivity_matrix(nodeNet)
    # locate the root node with largest radius
    # theoritically this should be the center of soma 
    _, rootNodeIndex = findmax(radii)
    seedNodeIndexes::Vector = [rootNodeIndex]
    while any(collectedFlagVec)
        if isempty(seedNodeIndexes)
            # this is not the first time, the main net is already built
            # find the uncollected node that is nearest to the branch list as seed 
            nearestNodeIndex, nearestBranchIndex, nearestNodeIndexInBranch = 
                        find_nearest_node_index(branchList, nodes, collectedFlagVec)
            seedNodeIndexes = [nearestNodeIndex]
        else 
            # if there still exist unvisitied edges
            while !isempty(seedIndexes)
                # grow the net from seed until there is branching point
                seedNodeIndex = pop!(seedIndexes)
                # grow a subnet
                subnet = BranchNet(seedNodeIndex, nodes, nodesConnectivityMatrix, collectedFlagVec)
                # merge the subnet to main net
                merge!(net, subnet, nearestBranchIndex, nearestNodeIndexInBranch)
            end 
        end 
    end 
    net     
end 

function BranchNet(seedNodeIndex::Integer, nodes::Vector{NTuple{4, Float32}}, 
             nodesConnectivityMatrix::SparseMatrixCSC, collectedFlagVec::Vector{Bool})
    branchList = Vector{Branch}()
    # the seed of branch should record both seed node index and the the parent branch index
    # the parent branch index of root branch is -1 
    branchSeedList = Vector{NTuple{Int, Int}}( [(seedNodeIndex, -1)] )
    # depth first search
    while !isempty(branchSeedList)
        seedNodeIndex, branchParentIndex = pop!(branchSeedList)
        # grow a branch from seed
        branch = Vector{NTuple{4, Float32}}()
        while true
            # construct this branch
            # push the seed node
            push!(branch, nodes[seedNodeIndex])
            # label this node index as collected
            collectedFlagVec[ seedNodeIndex ] = true 

            # find the connected nodes
            connectedNodeIndexList,_ = findnz(nodesConnectivityMatrix[:, seedNodeIndex])
            # exclude the collected nodes
            connectedNodeIndexList = connectedNodeIndexList[ !collectedFlagVec[connectedNodeIndexList] ] 

            if length(connectedNodeIndexList) == 0
                error("impossible, all node should connect to something!")
            elseif length(connectedNodeIndexList) == 1
                # belong to the same branch
                push!(branch, nodes[ connectedNodeIndexList[1] ])
                push!(seedNodeIndexList, connectedNodeIndexList[1])
            elseif length( connectedNodeIndexList ) > 1
                # multiple branching points, finish constructing this branch
                push!(branchList, branch)
                if branchParentIndex != -1
                    # this is not the root branch
                    connectivityMatrix[branchParentIndex, length(branchList)]
                end 
                # seed new branches
                for index in connectedNodeIndexList 
                    push!(branchSeedList, (index, length(branchList)))
                end 
                break
            end 
        end 
    end
    BranchNet(branchList, connectivityMatrix)
end 

"""
    get_num_branches(self::BranchNet)

get the number of branches  
"""
function get_num_branches(self::BranchNet) length(self.branchList) end 

"""
merge two nets at a specific branch location
"""
function Base.merge!(self::BranchNet, other::BranchNet, nearestBranchIndex::Integer, 
                     nearestNodeIndexInBranch::Integer)
    if isempty(self)
        self = other 
        return 
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
        self = BranchNet(self.branchList, mergedConnectivityMatrix)
        return
    end
end 

function Base.isempty(self::BranchNet)    isempty(self.branchList) end 


function NodeNets.NodeNet(net::BranchNet)
    
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

end # module
