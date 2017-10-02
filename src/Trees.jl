module Trees

using ..Skeletons

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
    start = convert(UInt32, start)
    
    mainTree = Tree( start, nodes, conn )
    # the remaining edges are disconnected with this tree 
    while !isempty( conn )

        for nodeIndex in keys(conn) 
            # find out the closest remaining point to grow a new tree

        end
        merge!(mainTree, newTree)
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

########## properties ######################
function get_root_branch(self::Tree) self.rootBranch end 
function get_radii(self::Tree) self.rootBranchRadii end 
function get_kind(self::Tree) self.kind end 
function get_children(self::Tree) self.children end 

############## Base functions ###################
# return the index of mainBranch and index of child
# it is a problem of how to iterate all the points in a tree!
# function Base.start( self::Tree )   1,0 end 

function Base.length( self::Tree )
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
"""
function Base.merge!(self::Tree,  point1::Union{Vector, Tuple}, 
                     other::Tree, point2::Union{Vector, Tuple})
    @assert length(point1) == length(point2) == 3
    error("unimplemented")
end 

############## measurements ####################
"""
compute the euclidean distance between a point and a tree 
"""
function euclidean(self::Tree, point::Union{Vector, Tuple})
    @assert length(point) == 3
    error("unimplemented")
end 

end # module
