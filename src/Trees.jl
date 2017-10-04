module Trees

using ..Skeletons
using ..SWCs

abstract AbstractTree 
export Tree 

const KIND = UInt8(0)

type Tree <: AbstractTree
    rootBranch          ::Array{UInt32, 2}
    rootBranchRadii     ::Vector{Float32} 
    children            ::Vector
    kind                ::UInt8
end 

function Tree( skeleton::Skeleton )
    nodes = Skeletons.get_nodes(skeleton)
    radii = Skeletons.get_radii(skeleton)
    # construct the connection dict to represent all the node connections
    conn = Skeletons.get_connection_dict(skeleton)

    # locate the root node with largest radius
    # theoritically this should be the soma center
    _, start = findmax( radii )
    root = convert(UInt32, start)
    
    mainTree = Tree( root, nodes, conn )
    # the remaining edges are disconnected with this tree 
    while !isempty( conn )
        # the minimum distance from a remaining node to the main tree
        minDistance = typemax(Float32)
        # the node index of nodes array, which is the closest to the mainTree
        minDistanceNodeIndex = UInt32(0)
        # the index of closest node. this index indicates to the main Tree
        closestMainTreeNodeIndex = UInt32(0)
        for nodeIndex in keys(conn) 
            # find out the closest remaining point to grow a new tree
            ramainingNode = nodes[nodeIndex]
            # distance from the mainTree
            distance, closestMainTreeNodeIndex = euclidean(mainTree, remainingNode)
            if distance < minDistance 
                minDistanceNodeIndex = nodeIndex 
            end 
        end
        newTree = Tree( minDistanceNodeIndex, nodes, conn )
        mainTreeNodes = get_nodes(mainTree)
        closestMainTreeNode = mainTreeNodes[ closestMainTreeNodeIndex, : ]
        # connect the closestMainTreeNode in main tree and the root node in the new tree.
        merge!(mainTree, closestMainTreeNode, newTree)
    end 
    return mainTree
end 

function Tree{T}( start::Integer, nodes::Array{T,2}, conn::Dict{T,Set{T}} )
    # construct the main branch until there is a branching node 
    rootBranch = Array{UInt32, 2}()
    rootBranchRadii = Vector{Float32}
    children = Vector{Tree}()
    while true
        if length(conn[start]) == 1
            # no branching point
            push!(rootBranch, nodes[start, :])
            push!(rootBranchRadii, radii[start])
            start = pop!(conn[start])
            delete!(conn, start)
        else 
            # multiple branches
            # push the branching point and start new tree from the branches
            push!(rootBranch, nodes[start, :])
            push!(rootBranchRadii, radii[start])
            for branchStartIndex in conn[start]
                push!(children, Tree( branchStartIndex, nodes, conn ))
                delete!(conn, branchStartIndex)
            end
            break
        end 
    end
    return Tree( rootBranch, rootBranchRadii, children, KIND ) 
end 

"""
transfer tree to swc 
"""
function SWC(tree::Tree; attributes = Dict{Symbol,Any}())
    rootBranch = Trees.get_root_branch( tree )

end 


########## properties ######################
function get_root_branch(self::Tree) self.rootBranch end 
function get_root_branch_radii(self::Tree) self.rootBranchRadii end 
function get_kind(self::Tree) self.kind end 
function get_children(self::Tree) self.children end 
"""
    get_nodes(self::Tree)

get all the nodes in the tree
"""
function get_nodes(self::Tree) 
    nodes = get_root_branch(self)
    for child in get_children(self)
        nodes = vcat(nodes, get_nodes(child))
    end 
    return nodes 
end 

############## Base functions ###################
# return the index of mainBranch and index of child
# it is a problem of how to iterate all the points in a tree!
# function Base.start( self::Tree )   1,0 end 

function get_total_length( self::Tree )
    l = Float32(0)
    # compute the root branch length
    for i in 1:length(self.rootBranch)-1
        l += norm(Vector{Float32}(self.rootBranch[i,  :]) .- 
                  Vector{Float32}(self.rootBranch[i+1,:]))
    end 
    # add the children tree length
    for child in self.children
        l += length( child )
    end 
    l
end

"""
merge two trees using the closest point pair
the root of other tree should connect to the closes point in the self tree  
"""
function Base.merge!(self::Tree,  closestNode::Vector, other::Tree )
    @assert length(point1) == length(point2) == 3
    isMerged = false
    rootBranch = get_root_branch(self)
    rootBranchRadii = get_root_branch_radii( self )
    kind = get_kind(self)
    children = get_children( self )

    for i in 1:size(rootBranch, 1)
        if all(rootBranch[i,:].== closestNode)
            # merge here
            newRootBranch = rootBranch[1:i, :]
            newRootBranchRadii = rootBranchRadii[1:i]
            if i < size(rootBranch, 1)
                newSubTree = Tree(rootBranch[i+1:end, :], rootBranchRadii[i+1:end], kind, children)
                newChildren = [newSubTree, other]
                self = Tree(newRootBranch, newRootBranchRadii, kind, newChildren)
                return true  
            else
                # merge to the branching point
                push!(self.children, other)
                return true 
            end 
        end 
    end 
    # if there is no matching node, keep merging
    for child in children 
        if merge!(child, closestNode, other)
            return true 
        end 
    end 
    return false 
end 

############## measurements ####################
"""
compute the euclidean distance between a point and a tree 
"""
function euclidean(self::Tree, point::Union{Vector, Tuple})
    @assert length(point) == 3
    nodes = get_nodes(self)
    euclidean(nodes, point)
end

function euclidean{T}(nodes::Array{T,2}, point::Union{Vector, Tuple})
    @assert length(point) == 3
    @assert size(nodes, 2) == 3
    ns = Array{Float32, 2}(nodes)
    minDistance = typemax(Float32)
    closestNodeIndex = 0
    p = Vector{Float32}(point)
    for i in 1:size(ns, 1)
        d = norm( ns[i,:] .- p )
        if d < minDistance 
            minDistance = d
            closestNodeIndex = i
        end 
    end 
    return distance, closestNodeIndex 
end 

end # module
