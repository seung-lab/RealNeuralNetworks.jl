module NBLASTs

using NearestNeighbors 
using LinearAlgebra 
using ProgressMeter 

using ..RealNeuralNetworks.Utils.VectorClouds 
using ..RealNeuralNetworks.NodeNets
using ..RealNeuralNetworks.SWCs 
using ..RealNeuralNetworks.Neurons 
using ..RealNeuralNetworks.Utils.Mathes
using ..RealNeuralNetworks.Utils.RangeIndexingArrays 
using ..RealNeuralNetworks.Neurons.Segments 


export nblast, nblast_allbyall


"""
    VectorCloud(neuron::SWC; k::Integer=20, class::Union{Nothing, UInt8}=nothing) 
"""
@inline function VectorCloud(neuron::SWC; k::Integer=20, 
                             class::Union{Nothing, UInt8}=nothing,
                             downscaleFactor::Int=1) 
    VectorCloud(NodeNet(neuron); k=k, class=class, downscaleFactor=downscaleFactor)
end 

"""
    VectorCloud(neuron::Neuron; k::Integer=20, class::Union{Nothing, UInt8}=nothing) 
"""
@inline function VectorCloud(neuron::Neuron{T}; k::Integer=20, 
                             class::Union{Nothing, UInt8}=nothing,
                             downscaleFactor::Int=1) where T
    VectorCloud(NodeNet(neuron); k=k, class=class, downscaleFactor=downscaleFactor)
end 

"""
    VectorCloud(neuron::NodeNet; k::Integer=20, class::Union{Nothing, UInt8}=nothing) 

Parameters: 
    neurons: input neuron 
    k: the number of neighboring nodes to use for principle vector computation
    class: the node type/class. The type/class was defined in SWC format.

Return: 
    ret::Matrix{Float32}: the vector cloud 
"""
function VectorCloud(neuron::NodeNet{T}; k::Integer=20, 
                    class::Union{Nothing, UInt8}=nothing,
                    downscaleFactor::Int=1) where T 
    nodeClassList = NodeNets.get_node_class_list(neuron) 

    # transform neuron to xyzmatrix
    N = NodeNets.get_node_num(neuron; class=class) 
    if N < k 
        @warn("the number of node $N is less than than the number of neighborhood points k=$k, return zero vector")
        return zeros(T, (0,0))
    end 

    xyzmatrix = Matrix{Float32}(undef, 3, N)

    if class == nothing 
        # all the nodes should be included 
        @inbounds for (i, node) in NodeNets.get_node_list(neuron) |> enumerate 
            xyzmatrix[:,i] = [node[1:3]...,]
        end 
    else 
        # only nodes match the class should be included 
        @assert isa(nodeClassList, Vector{UInt8})
        j = 0
        @inbounds for (i, node) in NodeNets.get_node_list(neuron) |> enumerate 
            if nodeClassList[i] == class 
                j += 1
                xyzmatrix[:, j] = [node[1:3]...,]
            end 
        end
    end 

    tree = KDTree(xyzmatrix; leafsize=k)

    # the first 3 rows will be X,Y,Z, and the last 3 rows will be X,Y,Z of direction vector
    ret = Matrix{Float32}(undef, 6, N)

    idxs, dists = knn(tree, xyzmatrix, k, false)

    data = Matrix{Float32}(undef, 3, k)
    @inbounds for (nodeIndex, indexList) in idxs |> enumerate
        data = xyzmatrix[:, indexList]
        PCs, eigenValueList = Mathes.pca(data)
        # use the first principle component as the direction vector
        ret[:, nodeIndex] = [xyzmatrix[:, nodeIndex]..., PCs[1]...,]
    end

    if downscaleFactor != 1
        ret[1:3,:] ./= T(1000)
    end
    ret 
end

function nblast(targetNeuron::Neuron{T}, queryNeuron::Neuron{T};
                downscaleFactor::Number=1,
                ria::RangeIndexingArray{T,2}=RangeIndexingArray{T}()) where T
    targetVectorCloud = VectorCloud(targetNeuron, downscaleFactor=downscaleFactor)
    queryVectorCloud = VectorCloud(queryNeuron, downscaleFactor=downscaleFactor)
    nblast(targetVectorCloud, queryVectorCloud; ria=ria)
end

function nblast(vectorCloudList::Vector{Matrix{T}}, 
                targetIndex::Integer, queryIndex::Integer;
                ria::RangeIndexingArray{T,2}=RangeIndexingArray{T}(),
                targetTree=VectorClouds.to_kd_tree(vectorCloud[targetIndex][1:3,:])) where {T}
    target = vectorCloudList[targetIndex]
    query = vectorCloudList[queryIndex]
    return nblast(target, query; ria=ria, targetTree=targetTree)
end

"""
    nblast( target::VectorCloud, query::VectorCloud;
            ria::RangeIndexingArray{T,2}=RangeIndexingArray{T,2}(),
            targetTree=KDTree(targetVectorCloud[1:3, :]))

measure the similarity of two neurons using NBLAST algorithm 

Reference:
```
Costa, Marta, et al. "NBLAST: rapid, sensitive comparison of neuronal structure and construction of neuron family databases." Neuron 91.2 (2016): 293-311.
```
"""
function nblast(target::Matrix{T}, query::Matrix{T}; 
                ria::RangeIndexingArray{T,2}=RangeIndexingArray{Float32}(), 
                targetTree::Union{Nothing, KDTree}=VectorClouds.to_kd_tree(target)) where T
   
    if isempty(target) || isempty(query) 
        # if one of them is empty, return zero 
        return zero(T)
    end 
    
    totalScore = zero(T)

    idxs, dists = knn(targetTree, query[1:3, :], 1, false)

    @inbounds for (i, nodeIndexList) in idxs |> enumerate
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
    normalize the raw NBLAST score by self score of query vector cloud 
"""
@inline function normalize_similarity_matrix!(similarityMatrix::Matrix{T}) where T
    @inbounds for i in 1:size(similarityMatrix, 1)
        similarityMatrix[:,i] ./= similarityMatrix[i,i]
    end
end

@inline function set_mean!(similarityMatrix::Matrix{T}) where T
    @inbounds for i in 1:size(similarityMatrix,1) 
        for j in i+1:size(similarityMatrix,2)
            similarityMatrix[i,j] = (similarityMatrix[i,j] + similarityMatrix[j,i])/T(2)
            similarityMatrix[j,i] = similarityMatrix[i,j]
        end 
    end 
end

"""
    nblast_allbyall(vectorCloudList::Vector{Union{Nothin, Matrix{T}}};
                    ria::RangeIndexingArray{TR,N}=RangeIndexingArray{Float32}(),
                    normalisation::Symbol=:raw) 

pairwise computation of similarity score and return a similarity matrix 

Parameters:
    ria: the precomputed score table, which as converted to RangeIndexingArray 
    vectorCloudList: a list of vector cloud for internal comparison 
    normalisation: normalisation method, could be one of raw|normalised|mean 
Return: 
    similarityMatrix::Matrix{TR}: the similarity matrix 
"""
function nblast_allbyall(vectorCloudList::Vector{Matrix{T}};
                        ria::RangeIndexingArray{T,2}=RangeIndexingArray{Float32}(),
                        treeList::Vector = map(VectorClouds.to_kd_tree, vectorCloudList)
                        ) where {T}
    num = length(vectorCloudList)
    similarityMatrix = Matrix{T}(undef, num, num)

    @inbounds @showprogress 1 "computing similarity matrix..." for targetIndex in 1:num 
    #Threads.@threads for targetIndex in 1:num 
        Threads.@threads for queryIndex in 1:num 
        #for queryIndex in 1:num 
            similarityMatrix[targetIndex, queryIndex] = nblast( 
                        vectorCloudList, targetIndex, queryIndex;
                        ria=ria, targetTree=treeList[targetIndex] )
        end 
    end 
   
    similarityMatrix
end

"""
    nblast_allbyall(neuronList::Vector{Neuron{T}}; 
                    semantic::Bool=false,
                    k::Int=20,
                    ria::RangeIndexingArray{TR}=RangeIndexingArray{Float32}(),
                    normalisation::Symbol=:raw) where {TR}
    Note that the neuron coordinate unit should be nm, it will be translated to micron internally.
    The semantic NBLAST will find the nearest vector pair with same semantic labeling. 
    An axonal vector in neuron A will find closest axonal vector in neuron B.
    A dendritic vector in neuron A will find closest dendritic vector in neuron B.
    It was implemented differently. We run two nblast computation with same semantic labeling 
    for raw score and sum the two scores up. Because the definition of raw score is a summation 
    of vector pair scores.

    Parameters:
        k: the number of nearest nodes used to compute the principle direction/vector
"""
function nblast_allbyall(neuronList::Vector{Neuron{T}};
                            semantic::Bool=false, 
                            k::Int=20,
                            ria::Union{Nothing, RangeIndexingArray{T,2}}=nothing,
                            downscaleFactor::Number=1000,
                            dendPositionOnly::Bool=true,
                            normalisation::Symbol=:raw) where T
    if ria == nothing 
        ria = RangeIndexingArray{T}()
    end
    if semantic
        # transforming to vector cloud list    
        vectorCloudListAxon = map(x->VectorCloud(x;class=Segments.AXON_CLASS,
                                    k=k, downscaleFactor=downscaleFactor), neuronList)
        similarityMatrixAxon = nblast_allbyall(vectorCloudListAxon; ria=ria)
        vectorCloudListDend = map(x->VectorCloud(x;class=Segments.DENDRITE_CLASS, 
                                    k=k, downscaleFactor=downscaleFactor), neuronList)
        if dendPositionOnly
            dendRia = RangeIndexingArrays.to_position_only(ria)
        else
            dendRia = ria
        end
        similarityMatrixDend = nblast_allbyall(vectorCloudListDend; ria=dendRia)
        similarityMatrix = similarityMatrixAxon .+ similarityMatrixDend
    else
        # transforming to vector cloud list    
        vectorCloudList = map(x->VectorCloud(x;class=nothing, 
                                    k=k, downscaleFactor=downscaleFactor), neuronList)
        similarityMatrix = nblast_allbyall(vectorCloudList; ria=ria)
    end
    
    if normalisation==:normalised || normalisation==:mean 
        normalize_similarity_matrix!(similarityMatrix)
    end 
    if normalisation==:mean
        set_mean!(similarityMatrix) 
    end

    return similarityMatrix
end


function small_to_big_nblast!(similarityMatrix, vectorCloudList, i, j, treeList, selfScoreList)
    vectorCloud1 = vectorCloudList[i]
    vectorCloud2 = vectorCloudList[j]
    # always use smaller neuron as query
    # normalize the raw score using query neuron
    if length(vectorCloud1) > length(vectorCloud2)
        score = nblast(vectorCloud1, vectorCloud2; ria=ria, targetTree=treeList[i]) / selfScoreList[j]
    else
        score = nblast(vectorCloud2, vectorCloud1; ria=ria, targetTree=treeList[j]) / selfScoreList[i]
    end
    similarityMatrix[i,j] = score
    similarityMatrix[j,i] = score
end

"""
    nblast_allbyall_small2big(neuronList::Vector{Neuron{T}}; 
            ria::Union{Nothing, RangeIndexingArray{T,2}}=nothing,
            treeList::Vector=map(VectorClouds.to_kd_tree, vectorCloudList),
            selfScoreList::Vector=map((v,t)->nblast(v,v; ria=ria, targetTree=t), 
                                                        vectorCloudList, treeList)) where T
"""
function nblast_allbyall_small2big(neuronList::Vector{Neuron{T}}; 
            ria::Union{Nothing, RangeIndexingArray{T,2}}=nothing,
            treeList::Vector=map(VectorClouds.to_kd_tree, vectorCloudList),
            selfScoreList::Vector=map((v,t)->nblast(v,v; ria=ria, targetTree=t), 
                                                        vectorCloudList, treeList)) where T
    vectorCloudList = map(x->VectorCloud(x; k=10, downscaleFactor=1000), neuronList);
    nblast_allbyall_small2big(vectorCloudList; ria=ria, treeList=treeList, 
                                                selfScoreList=selfScoreList)
end

"""
    nblast_allbyall_small2big(vectorCloudList::Vector{VectorCloud{T}}; 
            ria::Union{Nothing, RangeIndexingArray{T,2}}=nothing,
            treeList::Vector=map(VectorClouds.to_kd_tree, vectorCloudList),
            selfScoreList::Vector=map((v,t)->nblast(v,v; ria=ria, targetTree=t), 
                                                        vectorCloudList, treeList)) where T
"""
function nblast_allbyall_small2big(vectorCloudList::Vector{X}; 
            ria::Union{Nothing, RangeIndexingArray{T,2}}=nothing,
            treeList::Vector=map(VectorClouds.to_kd_tree, vectorCloudList),
            selfScoreList::Vector=map((v,t)->nblast(v,v; ria=ria, targetTree=t), 
                        vectorCloudList, treeList)) where {X<:Matrix{Float32}, T}
    N = length(vectorCloudList)
    similarityMatrix = ones(Float32, (N,N))
    # Threads.@threads for i in 1:N
    @showprogress for i in 1:N
        Threads.@threads for j in (i+1):N
    #     for j in (i+1):N
            small_to_big_nblast!(similarityMatrix, vectorCloudList, i, j, treeList, selfScoreList)
        end
    end
    similarityMatrix
end

end # end of module
