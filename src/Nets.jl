module Nets
include("Branches.jl")
using .Branches 
using ..TEASAR.Skeletons

export Net 
type Net 
    # x,y,z,r
    branchList ::Vector{NTuple{4, Float32}}
    connectivityMatrix ::SparseMatrixCSC{Bool, Int}
end 

function Net(skeleton::Skeleton)
    nodes = Skeletons.get_nodes(skeleton)
    radii = Skeletons.get_radii(skeleton)
    # locate the root node with largest radius
    # theoritically this should be the center of soma 
    _, start = findmax(radii)
end 

end # module
