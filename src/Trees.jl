module Trees

using ..Skeletons

abstract AbstractTree 
export Tree 

type Tree{T} <: AbstractTree
    mainBranch ::Array{T, 2}
    subTrees ::Vector 
end 

function Tree{T}( skeleton::Skeleton )
    
end 

end # module
