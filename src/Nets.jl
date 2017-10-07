module Nets
include("Branches.jl")
using .Branches 
using ..TEASAR.Skeletons

export Net 
type Net 
    # x,y,z,r
    branchList ::Vector{Branch}
    connectivityMatrix ::SparseMatrixCSC{Bool, Int}
end 

function Net(skeleton::Skeleton)
    numBranchPoint = Skeletons.get_num_branch_point(skeleton)
    # number of branches, assumes that this is acyclic graph
    numBranches = numBranchPoint + 2
    # initializae the net
    branchList = Vector{Branch}()
    connectivityMatrix = spzeros(Bool, numBranches, numBranches)
    net = Net(branchList, connectivityMatrix)

    # the properties from skeleton
    nodes = Skeletons.get_nodes(skeleton)
    radii = Skeletons.get_radii(skeleton)
    # flags labeling whether this node was collected to the net
    collectedFlagVec = zeros(Bool, length(nodes))
    # connectivity matrix of nodes in skeleton 
    nodesConnectivityMatrix  = Skeletons.get_connectivity_matrix(skeleton)
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
        end 
        # if there still exist unvisitied edges
        while !isempty(seedIndexes)
            # grow the net from seed until there is branching point
            seedNodeIndex = pop!(seedIndexes)
            # grow a subnet
            subnet = Net(seedNodeIndex, nodes, nodesConnectivityMatrix, collectedFlagVec)
            # merge the subnet to main net
            merge!(net, subnet, nearestBranchIndex, nearestNodeIndexInBranch)
        end 
    end 
    net     
end 

function Net(seedNodeIndex::Integer, nodes::Vector{NTuple{4, Float32}}, 
             nodesConnectivityMatrix::SparseMatrixCSC, collectedFlagVec::Vector{Bool})
    branchList = Vector{Branch}()
    seedIndexList = [seedNodeIndex]
    while !isempty(seedIndexList)
        seedIndex = pop!(seedIndexList)
        # grow a branch from seed
        branch = Vector{NTuple{4, Float32}}()
        while true 
            # push the seed node
            push!(branch, nodes[seedIndex])

            # find the connected nodes
            seedNodeIndexList,_ = findnz(nodesConnectivityMatrix[:, seedNodeIndex])
            # exclude the collected nodes
            seedNodeIndexList = 
            # 
            connectedNodes = 
            if length(seedNodeIndexes) == 0
                error("impossible, all node should connect to something!")
            elseif length(seedNodeIndexes) == 1
                # belong to the same 
                seedNodeIndex = seedIndexes[1]
                push!(branch, nodes[ seedNodeIndex ])
            elseif length(seedNodeIndexes) > 1
                # multiple 



end 

"""
    get_num_branches(self::Net)

get the number of branches  
"""
function get_num_branches(self::Net) length(self.branchList) end 

"""
merge two nets at a specific branch location
"""
function Base.merge!(self::Net, other::Net, nearestBranchIndex::Integer, 
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
        self = Net(self.branchList, mergedConnectivityMatrix)
        return
    end
end 

function Base.isempty(self::Net)    isempty(self.branchList) end 

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
