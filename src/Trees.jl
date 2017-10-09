module Trees

using ..NodeNets

abstract AbstractTree 
export Tree 

type Tree{T} <: AbstractTree
    mainBranch  ::Array{T, 2}
    radii       ::Vector 
    subTrees    ::Vector 
end 

function Tree{T}( nodeNet::NodeNet )
    
end 

function Base.length( self::Tree )

end 

end # module
