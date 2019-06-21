module VectorClouds

using LinearAlgebra
using NearestNeighbors
using ..Utils.Mathes 

const VectorCloud = Matrix

export VectorCloud

@inline function to_kd_tree(self::VectorCloud; leafsize::Integer=1)
    if isempty(self)
        return nothing 
    else 
        return KDTree(self[1:3, :], leafsize=leafsize)
    end 
end 


@inline function to_kd_tree(self::Nothing)
    return nothing
end

"""
    query(self::VectorCloud{T}, targetVectorCloud::VectorCloud; 
                targetTree::KDTree = to_kd_tree(targetVectorCloud))
query a VectorCloud to a target KDTree 
"""
function query(self::VectorCloud{T}, targetVectorCloud::VectorCloud{T}; 
               targetTree::KDTree=to_kd_tree(targetVectorCloud)) where T
    queryVectorCloud = self 

    idxs, dists = knn(targetTree, queryVectorCloud[1:3, :], 1, false)

    # distance matrix 
    ret = Matrix{T}(undef, 2, size(queryVectorCloud, 2))
    
    for (i, nodeIndexList) in idxs |> enumerate
        targetNodeIndex = nodeIndexList[1]
        # physical distance
        dist = dists[i][1]
        
        queryVector = queryVectorCloud[4:6, i]
        targetVector = targetVectorCloud[4:6, targetNodeIndex]
        # absolute dot product
        adp = abs( dot( queryVector, targetVector ) )
        
        ret[:, i] = [dist, adp]
    end
    ret
end 

end # end of module VectorClouds

