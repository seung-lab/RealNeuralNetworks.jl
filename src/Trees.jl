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
    
   return Tree( start, nodes, conn )
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
        else 
            # multiple branches
            # push the branching point and start new tree from the branches
            push!(rootBranch, nodes[start, :])
            push!(rootBranchRadii, radii[start])
            for branchStartIndex in conn[start]
                push!(children, Tree( branchStartIndex, nodes, conn ))
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
function Base.start( self::Tree )   self.rootBranch[1,:] end 

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


end # module
