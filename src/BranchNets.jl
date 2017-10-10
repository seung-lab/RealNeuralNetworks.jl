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
    # note that the connectivity matrix may not be a square matrix
    connectivityMatrix = sparse(parentBranchIndexList, childBranchIndexList, 
                                ones(Bool,length(childBranchIndexList)))
    BranchNet(branchList, connectivityMatrix)
end 

function BranchNet{T}( seg::Array{T,3}; obj_id::T = convert(T,1), 
                        expansion::NTuple{3,UInt32}=EXPANSION )
    nodeNet = NodeNet( seg; obj_id = obj_id, expansion = expansion )
    BranchNet( nodeNet )
end 

"""
    get_num_branches(self::BranchNet)

get the number of branches  
"""
function get_num_branches(self::BranchNet) length(self.branchList) end 
function get_branch_list(self::BranchNet) self.branchList end 
function get_connectivity_matrix(self::BranchNet) self.connectivityMatrix end 
function get_num_node_edges(self::BranchNet) error("unimplemented") end 

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
        # just stitch the new branches and rebuild the connection matrix
        mergedBranchList = vcat( branchList1, branchList2 )
        # total number of branches
        total_num_branches = num_branches1 + num_branches2 
        mergedConnectivityMatrix = spzeros(Bool, total_num_branches, total_num_branches)
        mergedConnectivityMatrix[1:size(self.connectivityMatrix,1), 
                                 1:size(self.connectivityMatrix,2)] = 
                                                                self.connectivityMatrix 

        mergedConnectivityMatrix[
                num_branches1+1 : num_branches1+1+size(other.connectivityMatrix,1)-1, 
                num_branches1+1 : num_branches1+1+size(other.connectivityMatrix,2)-1] = 
                                                                other.connectivityMatrix

        # connect the two nets, the row index is the parent and the column is the child
        mergedConnectivityMatrix[nearestBranchIndex, num_branches1+1] = true
        return BranchNet(mergedBranchList, mergedConnectivityMatrix)
    else 
        # need to break the nearest branch and rebuild connectivity matrix
        total_num_branches = num_branches1 + num_branches2 + 1
        mergedBranchList = branchList1 
        mergedConnectivityMatrix = spzeros(Bool, total_num_branches, total_num_branches)
        
        # need to break the branch and then stitch the new branches
        branchPart1, branchPart2 = split(self.branchList[nearestBranchIndex], 
                                                nearestNodeIndexInBranch)
        mergedBranchList[nearestBranchIndex] = branchPart1 
        mergedConnectivityMatrix[1:size(self.connectivityMatrix,1), 
                                 1:size(self.connectivityMatrix,2)] = 
                                                        self.connectivityMatrix
        # reconnect the breaked two branches
        push!(mergedBranchList, branchPart2)
        mergedConnectivityMatrix[nearestBranchIndex, num_branches1+1] = true 

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

########################## type convertion ####################
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


"""
    find_nearest_node_index(branchList::Vector{Branch}, nodes::Vector{NTuple{4,Float32}})

find the uncollected node which is nearest to the branch list 
"""
function find_nearest_node_index(branchList::Vector{Branch}, 
                                 nodes::Vector{NTuple{4,Float32}},
                                 collectedFlagVec ::Vector{Bool})
    ret = (0,0,0)
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
    # the branchIndex should be inside the branchList
    @assert ret[2] <= length(branchList)
    @assert ret!=(0,0,0) 
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
