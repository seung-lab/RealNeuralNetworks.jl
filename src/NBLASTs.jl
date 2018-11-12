module NBLASTs

using ..RealNeuralNetworks.Utils.VectorClouds 
using ..RealNeuralNetworks.NodeNets
using ..RealNeuralNetworks.SWCs 
using ..RealNeuralNetworks.Neurons 
using ..RealNeuralNetworks.Utils.Mathes
using ..RealNeuralNetworks.Utils.RangeIndexingArrays 

using NearestNeighbors 
using LinearAlgebra 
using ProgressMeter 

export nblast, nblast_allbyall

@inline function VectorCloud(neuron::SWC; k::Integer=20)
    VectorCloud(NodeNet(neuron); k=k)
end 

@inline function VectorCloud(neuron::Neuron; k::Integer=20)
    VectorCloud(NodeNet(neuron); k=k)
end 

"""
    VectorCloud(neuron::NodeNet; k::Integer=20) 
Parameters: 
    neurons: input neuron 
    k::Integer: the number of neighboring nodes to use for principle vector computation

Return: 
    ret::Matrix{Float32}: the vector cloud 
"""
function VectorCloud(neuron::NodeNet; k::Integer=20) 
    # transform neuron to xyzmatrix
    xyzmatrix = Matrix{Float32}(undef, 3, NodeNets.get_num_nodes(neuron))

    for (j, node) in NodeNets.get_node_list(neuron) |> enumerate
        xyzmatrix[:, j] = [node[1:3]...,]
    end
    tree = KDTree(xyzmatrix; leafsize=k)

    # the first 3 column will be X,Y,Z, and the last 3 columns will be X,Y,Z of direction vector
    ret = Matrix{Float32}(undef, 6, NodeNets.get_num_nodes(neuron))

    idxs, dists = knn(tree, xyzmatrix, k, false)

    data = Matrix{Float32}(undef, 3, k)
    for (nodeIndex, indexList) in idxs |> enumerate
        data = xyzmatrix[:, indexList]
        PCs, eigenValueList = Mathes.pca(data)
        # use the first principle component as the direction vector
        ret[:, nodeIndex] = [xyzmatrix[:, nodeIndex]..., PCs[1]...,]
    end
    ret 
end 

"""
    nblast(ria::RangeIndexingArray, target::VectorCloud, query::VectorCloud; 
                                                targetTree=KDTree(targetVectorCloud[1:3, :]))

measure the similarity of two neurons using NBLAST algorithm 
Reference:
Costa, Marta, et al. "NBLAST: rapid, sensitive comparison of neuronal structure and construction of neuron family databases." Neuron 91.2 (2016): 293-311.
"""
function nblast(ria::RangeIndexingArray{T,N}, 
                target::Matrix{T}, query::Matrix{T}; 
                targetTree::KDTree=VectorClouds.to_kd_tree(target)) where {T,N}
    
    totalScore = zero(Float32)

    idxs, dists = knn(targetTree, query[1:3, :], 1, false)

    for (i, nodeIndexList) in idxs |> enumerate
        targetNodeIndex = nodeIndexList[1]
        # physical distance
        dist = dists[i][1]

        queryVector = query[4:6, i]
        targetVector = target[4:6, targetNodeIndex]
        # absolute dot product
        adp = abs( dot( queryVector, targetVector ) )

        totalScore += ria[dist, adp]
    end
    totalScore
end

"""
    nblast_allbyall(ria::RangeIndexingArray{TR, N}, vectorCloudList::Vector{VectorCloud{T}};
                         normalisation::Symbol=:raw) 

pairwise computation of similarity score and return a similarity matrix 

Parameters:
    ria: the precomputed score table, which as converted to RangeIndexingArray 
    vectorCloudList: a list of vector cloud for internal comparison 
    normalisation: normalisation method, could be one of raw|normalised|mean 
Return: 
    similarityMatrix::Matrix{TR}: the similarity matrix 
"""
function nblast_allbyall(ria::RangeIndexingArray{TR, N}, 
                         vectorCloudList::Vector{Matrix{T}};
                         normalisation::Symbol=:raw) where {TR, N, T}
    num = length(vectorCloudList)
    similarityMatrix = Matrix{TR}(undef, num, num)

    treeList = map(VectorClouds.to_kd_tree, vectorCloudList) 

    @showprogress 1 "computing similarity matrix..." for targetIndex in 1:num 
        Threads.@threads for queryIndex in 1:num 
            similarityMatrix[targetIndex, queryIndex] = nblast(ria, 
                    vectorCloudList[targetIndex], vectorCloudList[queryIndex]; 
                    targetTree=treeList[targetIndex] )
        end 
    end 
   
    if normalisation==:normalised || normalisation==:mean 
        for i in 1:num
            similarityMatrix[:,i] ./= similarityMatrix[i,i]
        end
    end 
    if normalisation==:mean 
        for i in 1:num 
            for j in i+1:num 
                similarityMatrix[i,j] = (similarityMatrix[i,j] + similarityMatrix[j,i])/Float32(2)
                similarityMatrix[j,i] = similarityMatrix[i,j]
            end 
        end 
    end
    similarityMatrix
end 

end # end of module
