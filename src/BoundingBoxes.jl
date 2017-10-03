module BoundingBoxes 
type BoundingBox{T}
    minCorner ::Vector{T}
    maxCorner ::Vector{T} 

    function (::Type{BoundingBox}){T}( minCorner::Vector{T}, maxCorner::Vector{T})
        new{T}(minCorner, maxCorner)
    end 
end 




end # module
