using SparseArrays 

export NodeNet 
mutable struct NodeNet{T} 
    # the classes following the definition of swc 
    # 0 - undefined
    # 1 - soma
    # 2 - axon
    # 3 - (basal) dendrite
    # 4 - apical dendrite
    # 5 - fork point
    # 6 - end point
    # 7 - custom
    classes             :: Vector{UInt8}
    # each column: x,y,z,r
    # the size should be (4, nodeNum) 
    nodeArray           :: Matrix{T}
    # connectivity matrix to represent edges
    # conn[2,3]=true means node 2's parent is node 3
    # we have not using an array to store family relationship 
    # because we are assuming arbitrary trees with multiple childrent!
    connectivityMatrix  :: SparseMatrixCSC{Bool,UInt32}
end 
