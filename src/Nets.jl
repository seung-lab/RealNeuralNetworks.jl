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
    # initializae the net
    branchList = Vector{Branch}()
    connectivityMatrix = 
    # the properties from skeleton
    nodes = Skeletons.get_nodes(skeleton)
    radii = Skeletons.get_radii(skeleton)
    conn  = Skeletons.get_connectivity_matrix(skeleton)
    # locate the root node with largest radius
    # theoritically this should be the center of soma 
    _, root = findmax(radii)
    seeds::Vector = [root]
    while countnz(conn) > 0
        # if there still exist unvisitied edges
        while !isempty(seeds)
            # grow the net from seed until there is branching point
            seed = pop!(seeds)
            
end 

end # module
