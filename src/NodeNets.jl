export load_swc, save_swc, get_neuroglancer_precomputed, save_nodenet_bin, load_nodenet_bin
export stretch_coordinates!, add_offset!
export set_radius!
export get_path_length, downsample!, find_closest_node_id, get_radii

include("TEASAR.jl")

import LinearAlgebra: norm
import DelimitedFiles: readdlm, writedlm
using NearestNeighbors


const OFFSET = (zero(UInt32), zero(UInt32), zero(UInt32))
# rescale the skeleton
const EXPANSION = (one(UInt32), one(UInt32), one(UInt32))


@inline function Matrix{T}(self::NodeNet) where T
    parents = connectivity_matrix_to_parents(self.connectivityMatrix)
    hcat(self.classes, self.nodeArray, parents)
end

function NodeNet(nodeArray::Matrix{T}, 
                    connectivityMatrix::SparseMatrixCSC{Bool, UInt32}) where T
    @assert size(nodeArray, 1) == 4
    classes = zeros(UInt8, size(nodeArray, 2))
    dropzeros!(connectivityMatrix)
    NodeNet(nodeArray, classes, connectivityMatrix)
end

@inline function NodeNet(nodeList::Vector{NTuple{4,T}}, 
                         connectivityMatrix::SparseMatrixCSC{Bool, UInt32}) where {T}
    nodeArray = node_list_to_array(nodeList)
    NodeNet{T}(nodeArray, connectivityMatrix)
end 

"""
    NodeNet( seg, obj_id; penalty_fn=alexs_penalty)
Perform the teasar algorithm on the passed binary array.
"""
function NodeNet( seg::Array{T,3}; 
                     obj_id::T = convert(T,1), 
                     expansion::NTuple{3, UInt32} = EXPANSION,
                     penalty_fn::Function = alexs_penalty ) where T
    # note that the object voxels are false and non-object voxels are true!
    # bin_im = DBFs.create_binary_image( seg, obj_id ) 
    points = PointArrays.from_seg(seg; obj_id=obj_id)
    teasar(points; expansion=expansion, penalty_fn=penalty_fn) 
end 

"""
    NodeNet(bin_im)
Parameters:
    bin_im: binary mask. the object voxel should be false, non-object voxel should be true
Return:
    nodeNet object
"""
function NodeNet(bin_im::Union{BitArray, Array{Bool,3}}; 
                 offset::NTuple{3, UInt32} = OFFSET,
                 expansion::NTuple{3, UInt32} = EXPANSION,
                 penalty_fn::Function = alexs_penalty)
        # transform segmentation to points
    points = PointArrays.from_binary_image(bin_im)
    
    println("computing DBF");
    # boundary_point_indexes = PointArrays.get_boundary_point_indexes(points, seg; obj_id=obj_id)
    #@time DBF = DBFs.compute_DBF( points, boundary_point_indexes );
    @time DBF = DBFs.compute_DBF(points)
    # @time dbf = DBFs.compute_DBF(points, bin_im)

    PointArrays.add_offset!(points, offset)
    teasar(points; dbf=dbf, penalty_fn=penalty_fn, expansion = expansion)
end 


"""
Note that the root node id is 0 rather than -1
"""
function parents_to_connectivity_matrix(parents::Vector{UInt32})
    nodeNum = length(parents)

    childNodeIdxList = Vector{UInt32}()
    sizehint!(childNodeIdxList, nodeNum)
    parentNodeIdxList = Vector{UInt32}()
    sizehint!(parentNodeIdxList, nodeNum)

    @inbounds for childNodeIdx in 1:nodeNum
        parentNodeIdx = parents[childNodeIdx]
        if parentNodeIdx > zero(UInt32)
            push!(childNodeIdxList, childNodeIdx)
            push!(parentNodeIdxList, parentNodeIdx)
        end
    end

    connectivityMatrix = sparse(childNodeIdxList, parentNodeIdxList, 
                                                true, nodeNum, nodeNum)
    connectivityMatrix
end

function connectivity_matrix_to_parents(connectivityMatrix::SparseMatrixCSC{Bool, UInt32})
    nodeNum = size(connectivityMatrix, 1)
    parents = zeros(UInt32, nodeNum)
    
    childNodeIdxes, parentNodeIdxes, _ = findnz(connectivityMatrix)
    @inbounds for (childNodeIdx, parentNodeIdx) in zip(childNodeIdxes, parentNodeIdxes)
        parents[ childNodeIdx ] = parentNodeIdx
    end  
    parents
end

##################### getters ###############################
@inline function get_node_array(self::NodeNet) self.nodeArray end 

@inline function get_node_list(self::NodeNet{T}) where T 
    nodeArray = get_node_array(self)
    nodeNum = size(nodeArray, 2)
    nodeList = Vector{NTuple{4,T}}()
    sizehint!(nodeList, nodeNum)
    @inbounds for i in 1:nodeNum
        node = tuple(view(nodeArray, :, i)...)
        push!(nodeList, node)
    end
    nodeList
end

@inline function get_connectivity_matrix(self::NodeNet) self.connectivityMatrix end
@inline function get_classes(self::NodeNet) self.classes end 
@inline function get_radii(self::NodeNet) @view self.nodeArray[4, :] end 
function get_node_num(self::NodeNet; class::Union{Nothing, UInt8}=nothing)
    if class == nothing 
        return length(self.classes) 
    else 
        return count(self.classes .== class)
    end 
end

@inline function get_parents(self::NodeNet)
    connectivityMatrix = get_connectivity_matrix(self)
    connectivity_matrix_to_parents(connectivityMatrix)
end

@inline function get_root_node_index_list(self::NodeNet)
    parents = get_parents(self)
    findall(parents.==zero(UInt32))
end

function get_children_node_index_list(self::NodeNet, nodeIdx::Integer)
    connectivityMatrix = get_connectivity_matrix(self)
    childrenNodeIdxes,_ = findnz(connectivityMatrix[:, nodeIdx])
    childrenNodeIdxes
end

function get_parent_node_index(self::NodeNet, nodeIdx::Integer)
    connectivityMatrix = get_connectivity_matrix(self)
    parentNodeIdxes,_ = findnz(connectivityMatrix[:, nodeIdx])
    if length(parentNodeIdxes) == 1
        return parentNodeIdxes[1]
    elseif isempty(parentNodeIdxes)
        return zero(UInt32)
    else
        error("you should not have multiple parents: ", parentNodeIdxes)
    end
end

"""
    node_list_to_array(nodeList::Vector{NTuple{4,T}}) where T

transform a list of nodes to array. The size of array is (4, N).
The first axis is the x,y,z,r, the second axis is the nodes.
"""
function node_list_to_array(nodeList::Vector{NTuple{4,T}}) where T
    nodeNum = length(nodeList)
    ret = Array{T,2}(undef, 4, nodeNum)
    @inbounds for (i,node) in enumerate(nodeList)
        ret[:, i] = [node...]
    end
    ret
end

""" 
the connectivity matrix is symmetric, so the connection is undirected
"""
@inline function get_edge_num(self::NodeNet) nnz(self.connectivityMatrix) end

function get_terminal_node_id_list(self::NodeNet)
    terminalNodeIdList = Vector{Int}()
    connectivityMatrix = get_connectivity_matrix(self)
    dropzeros!(connectivityMatrix)
    for i in 1:get_node_num(self)
        if nnz(connectivityMatrix[:, i]) == 0
            push!(terminalNodeIdList, i)
        end 
    end
    return terminalNodeIdList
end 

function get_branching_node_id_list(self::NodeNet)
    branchingNodeIdList = Vector{Int}()
    connectivityMatrix = get_connectivity_matrix(self)
    dropzeros!(connectivityMatrix)
    for i in 1:get_node_num(self)
        if nnz(connectivityMatrix[:, i]) > 1
            push!(branchingNodeIdList, i)
        end 
    end
    return branchingNodeIdList
end 

#################### Setters ############################################
function set_radius!(self::NodeNet, radius::Float32)
    nodeArray = get_node_array(self)
    nodeArray[4, :] .= radius
end

#################### Base functions ######################################
function Base.isequal(self::NodeNet{T}, other::NodeNet{T}) where T
    parents1 = get_parents(self)
    parents2 = get_parents(other) 
    all(self.classes .== other.classes) && 
    all(self.nodeArray .== other.nodeArray) && 
    all(parents1 .== parents2)
end

@inline function Base.:(==)(self::NodeNet{T}, other::NodeNet{T}) where T
    isequal(self, other)
end

@inline function Base.length(self::NodeNet)
    length(self.classes)
end 

"""
    Base.iterate(self::NodeNet, state::NTuple{2,Int}=(1,1))

depth first search of tree. 
We assume that nodeNet is a forest of trees.
The first element of state is tree index.
The second element of state is node index in forrest (not in tree!).
"""
function Base.iterate(self::NodeNet, state::NTuple{2,Int}=(1,1); root)
    treeIdx, nodeIdx = state
    error("unimplemented")
end

############################### Base functions #######################################

"""
    Base.getindex(self::NodeNet, idx::Integer)

return a node of this neuron. The node is a view without memory copy.
The vector is [x,y,z,r] 
"""
@inline function Base.getindex(self::NodeNet, idx::Integer)
    nodeArray = get_node_array(self)
    @view nodeArray[:, idx]
end

"""
    Base.getindex(self::NodeNet, idx::Union{UnitRange, Vector{Int}})

Create a new NodeNet with a subset of current nodeNet
"""
function Base.getindex(self::NodeNet, selectedNodeIdxes::Union{UnitRange, Vector{Int}})
    error("this is not working correctly now.")
    nodeNum = length(self)
    newClasses = get_classes(self)[selectedNodeIdxes]
    newNodeArray = get_node_array(self)[:, selectedNodeIdxes]

    newNodeNum = length( selectedNodeIdxes )
    
    newParents = zeros(UInt32, newNodeNum)
    if isa(selectedNodeIdxes, Vector)
        selectedNodeIdxes = Set( selectedNodeIdxes )
    end

    for (newIdx, oldIdx) in enumerate( selectedNodeIdxes )
        nodeIdx = oldIdx
        parentNodeIdx = get_parent_node_index(self, nodeIdx)

        # find the parent/grandpar node and connect them
        while !(parentNodeIdx in idx)
            nodeIdx = parentNodeIdx
            parentNodeIdx = get_parent_node_index(self, nodeIdx)
        end
        newParents[newIdx] = parentNodeIdx
    end
    
    # map the parents index from old to new
    update_parents!(newParents, idx)

    newConnectivityMatrix = parents_to_connectivity_matrix(newParents)
    NodeNet(newClasses, newNodeArray, newConnectivityMatrix)
end 

function Base.UnitRange(self::NodeNet)
    nodeArray = get_node_array(self)
    minCoordinates = [typemax(UInt32), typemax(UInt32), typemax(UInt32)]
    maxCoordinates = [zero(UInt32), zero(UInt32), zero(UInt32)]
    @inbounds for i in 1:length(self)
        node = @view nodeArray[1:3, i]
        minCoordinates = map(min, minCoordinates, node)
        maxCoordinates = map(max, maxCoordinates, node)
    end 
    return [minCoordinates[1]:maxCoordinates[1], 
            minCoordinates[2]:maxCoordinates[2], 
            minCoordinates[3]:maxCoordinates[3]]
end 

"""
    find_closest_node_id(self::NodeNet{T}, point::NTuple{3,T}) where T

look for the id of the closest node
"""
@inline function find_closest_node_id(self::NodeNet{T}, point::NTuple{N,T}) where {N,T}
    find_closest_node_id(self, [point[1:3]...])
end

function find_closest_node_id(self::NodeNet{T}, point::Vector{T}) where T
    nodeArray = get_node_array(self)
    kdtree = KDTree(nodeArray[1:3, :]; leafsize=10)
    idxs, _ = knn(kdtree, point, 1)
    return idxs[1]
end

"""
    get_path_length( self::NodeNet )
accumulate all the euclidean distance of edges 
"""
function get_path_length( self::NodeNet{T} ) where T
    nodeArray = get_node_array(self)
    connectivityMatrix = get_connectivity_matrix(self)
    childNodeIdxes, parentNodeIdxes, _ = findnz(connectivityMatrix)

    pathLength = zero(T)
    @inbounds for (childNodeIdx, parentNodeIdx) in zip(childNodeIdxes, parentNodeIdxes)
        childNode = view(nodeArray, 1:3, childNodeIdx)
        parentNode = view(nodeArray, 1:3, parentNodeIdx)
        pathLength += norm(childNode .- parentNode)
    end
    pathLength
end 

##################### transformation ##########################
"""
get binary buffer formatted as neuroglancer nodeNet.

# Binary format
    UInt32: number of vertex
    UInt32: number of edges
    Array{Float32,2}: Nx3 array, xyz coordinates of vertex
    Array{UInt32,2}: Mx2 arrray, node index pair of edges
reference: 
https://github.com/seung-lab/neuroglancer/wiki/Skeletons

TO-DO:
Will have saved the radius as attributes, we can try to read that.
"""
function get_neuroglancer_precomputed(self::NodeNet)
    nodeNum = get_node_num(self)
    edgeNum = get_edge_num(self)
    # total number of bytes
    byteNum = 4 + 4 + 4*3*nodeNum + 4*2*edgeNum
    io = IOBuffer( read=false, write=true, maxsize=byteNum )

    # write the number of vertex, and edges
    write(io, UInt32(nodeNum))
    write(io, UInt32(edgeNum))
    
    # write the node coordinates
    nodeArray = get_node_array(self)
    write(io, nodeArray[1:3, :])
    
    # write the edges
    connectivityMatrix = get_connectivity_matrix(self)
    childNodeIdxes, parentNodeIdxes, _ = findnz(connectivityMatrix)
    edges = Matrix{UInt32}(undef, 2, edgeNum)
    # neuroglancer index start from 0
    edges[1, :] = childNodeIdxes .- one(UInt32)
    edges[2, :] = parentNodeIdxes .- one(UInt32)
    write(io, edges)
    
    data = take!(io)
    close(io)
    return data
end 

function deserialize(data::Vector{UInt8})
    # a pointObj is 21 byte
    @assert mod(length(data), 21) == 0 "the binary file do not match the byte layout of pointObj."
    nodeNum = div(length(data), 21)
    classes = Vector{UInt8}(undef, nodeNum)
    nodeArray = Matrix{Float32}(undef, 4, nodeNum)
    parents = zeros(UInt32, nodeNum) 

    @inbounds for i in 1:nodeNum
        nodeData = view(data, (i-1)*21+1 : i*21)
        classes[i] = nodeData[1]
        nodeArray[:, i] = reinterpret(Float32, nodeData[2:17])
        parents[i] = reinterpret(UInt32, nodeData[18:21])[1] 
    end 
    connectivityMatrix = parents_to_connectivity_matrix(parents)    
    
    NodeNet(classes, nodeArray, connectivityMatrix)
end

"""
    load_nodenet_bin( fileName::AbstractString )
"""
function load_nodenet_bin( fileName::AbstractString )
    @assert endswith(fileName, ".nodenet.bin")
    read( fileName ) |> deserialize
end 

function serialize(self::NodeNet)
    classes = get_classes(self)
    nodeArray = get_node_array(self)
    connectivityMatrix = get_connectivity_matrix(self)
    parents = connectivity_matrix_to_parents(connectivityMatrix)

    nodeNum = length(classes)
    byteNum = nodeNum * 21
    io = IOBuffer( read=false, write=true, maxsize=byteNum )
    for i in 1:nodeNum
        write(io, classes[i])
        write(io, nodeArray[:, i])
        write(io, parents[i])
    end
    
    data = take!(io)  
    close(io)
    @assert length(data) == byteNum
    data 
end 

"""
    save_swc_bin( self::NodeNet, fileName::AbstractString )
represent swc file as binary file. the data structure is the same with swc.
"""
function save_nodenet_bin( self::NodeNet, fileName::AbstractString )
    @assert endswith(fileName, ".nodenet.bin")
    data = serialize(self)
    write(fileName, data)
end 

@inline function load_swc(fileName::AbstractString)
    data = readdlm(fileName, ' ', Float32, '\n', comments=true, comment_char='#')
    @assert size(data, 2) == 7
    # the node ID should in order, so we can ignore this redundent information
    data = data[sortperm(data[:, 1]), :]
    # the index should be 1 based and continous
    @assert data[end, 1] == size(data, 1)

    classes = UInt8.(data[:, 2])
    nodeArray = data[:, 3:6]'|>Matrix

    # the root node id should be zero rather than -1 in our data structure 
    parents = @view data[:, 7]
    parents[parents.<zero(Float32)] .= zero(Float32)
    parents = UInt32.(parents)
    connectivityMatrix = parents_to_connectivity_matrix(parents)
    NodeNet(classes, nodeArray, connectivityMatrix)
end

"""
current implementation truncate the value to digits 3!
The integration transformation will loos some precision!
Currently, it is ok because our unit is nm and the resolution is high enough.
If it becomes a problem, we can use list of tuple, and the types are mixed in the tuple.

```julia
x = [1; 2; 3; 4];
y = [5.2; 6.3; 7.5; 8.7];
open("delim_file.txt", "w") do io
    writedlm(io, map(identity, zip(x,y)), ' ')
end
```
"""
function save_swc(self::NodeNet, file_name::AbstractString; truncDigits::Int=3)
    classes = get_classes(self)
    nodeArray = get_node_array(self)
    parents = get_parents(self)

    truncedNodeArray = trunc.(nodeArray; digits=truncDigits)
    xs = view(truncedNodeArray, 1, :)
    ys = view(truncedNodeArray, 2, :)
    zs = view(truncedNodeArray, 3, :)
    rs = view(truncedNodeArray, 4, :) 

    nodeNum = length(classes)
    data = zip(1:nodeNum, classes, xs, ys, zs, rs, parents)
    writedlm(file_name, data, ' ', dims=(nodeNum, 7))
end

##################### manipulate ############################

@inline function squred_distance(node1, node2)
    @fastmath   (node1[1] - node2[1]) * (node1[1] - node2[1]) + 
                (node1[2] - node2[2]) * (node1[2] - node2[2]) +
                (node1[3] - node2[3]) * (node1[3] - node2[3])
end

"""
    downsample!(self::NodeNet{T}; step::T = Float32(1000))

get the mean coordinate and radius of neighboring nodes in a segment 
replace the nearby nodes.
The default is 1 micron.
"""
function downsample!(self::NodeNet{T}; step::T = Float32(1000)) where T
    nodeNum = length(self)
    nodeArray = get_node_array(self)
    
    # this is used as a threshold, so we do not need to compute the sqrt
    stepSquared = step * step

    # use the root nodes as seeds
    seedNodeIdxes = UInt32.(get_root_node_index_list(self))
    seedNodeParentIdexes = zeros(UInt32, length(seedNodeIdxes))

    # we select some nodes out as new nodeNet 
    # this selection should include all the seed node indexes
    selectedNodeIdxes = UInt32[]
    # we need to update the parents when we add new seletected nodes.
    newParents = UInt32[]

    @fastmath @inbounds while !isempty(seedNodeIdxes)
        seedNodeIdx = pop!(seedNodeIdxes)
        parentNodeIdx = pop!(seedNodeParentIdexes)
        push!(selectedNodeIdxes, seedNodeIdx)
        push!(newParents, parentNodeIdx)

        childrenNodeIdxes = get_children_node_index_list(self, seedNodeIdx)
        childrenNum = length(childrenNodeIdxes)

        # measure the walking distance 
        walkDistance = zero(T)
        startNodeIdx = seedNodeIdx
        startNode = @view nodeArray[:, startNodeIdx]
        # the parent node index of current search

        # walk through a segment
        while childrenNum == 1
            # this node inside a segment 
            nodeIdx = childrenNodeIdxes[1]
            node = @view nodeArray[:, nodeIdx]
            childrenNodeIdxes = get_children_node_index_list(self, nodeIdx)
            childrenNum = length(childrenNodeIdxes)
            d2 = squred_distance(startNode, node)
            if d2 < stepSquared
                # continue search                    
                # current node becomes parent and will give birth of children
                parentNodeIdx = nodeIdx
            else
                # have enough walking distance, will include this node in new nodeNet
                push!(selectedNodeIdxes, nodeIdx)
                push!(newParents, startNodeIdx)

                # now we walk start from here
                startNodeIdx = nodeIdx
                startNode = @view nodeArray[:, startNodeIdx]
               
                # adjust the coordinate and radius by mean of nearest nodes.
                parentNode = @view nodeArray[:, parentNodeIdx]
                if childrenNum == 1
                    childNode = @view nodeArray[:, childrenNodeIdxes[1]]
                    nodeArray[:, nodeIdx] = (parentNode .+ node .+ childNode) ./ T(3)
                else
                    # this node is the end of segment
                    nodeArray[:, nodeIdx] = (parentNode .+ node) ./ T(2)
                end
            end
        end
        # add all children nodes as sees.
        # if reaching the terminal and there is no children
        # nothing will be added 
        append!(seedNodeIdxes, childrenNodeIdxes)
        for idx in childrenNodeIdxes
            push!(seedNodeParentIdexes, parentNodeIdx)
        end
    end
    
    update_parents!(newParents, selectedNodeIdxes)

    # pickout the selected nodes as new nodeNet
    newClasses = get_classes(self)[selectedNodeIdxes]
    newNodeArray = nodeArray[:, selectedNodeIdxes]
    newConnectivityMatrix = parents_to_connectivity_matrix(newParents)
    NodeNet(newClasses, newNodeArray, newConnectivityMatrix)
end

@inline function add_offset!(self::NodeNet{T}, offset::NTuple{3,T} ) where T
    add_offset!(self, [offset...])
end

function add_offset!(self::NodeNet{T}, offset::Vector{T} ) where T
    @assert length(offset) == 3
    nodeArray = get_node_array(self)
    @inbounds for i in size(nodeArray, 2)
        nodeArray[1:3, i] .+= offset
    end
end

@inline function stretch_coordinates!(self::NodeNet, mip::Real)
    expansion = [2^(mip-1), 2^(mip-1), 1]
    stretch_coordinates!(self, expansion)
end 

function stretch_coordinates!(self::NodeNet{T}, expansion::Union{Vector, Tuple}) where T
    @assert length(expansion) == 3
    radiusExpansion = prod(expansion)^(1/3)
    expansion = T.([expansion..., radiusExpansion])

    nodeArray = get_node_array(self)
    @inbounds for i in 1:size(nodeArray, 2)
        nodeArray[:, i] .*= expansion
    end
end

################################ utilities #################################################
function update_parents!(parents::Vector{UInt32}, selectedNodeIdxes::Union{UnitRange, Vector})
    # this is a map connecting old and new index
    oldIdx2newIdx = zeros(UInt32, maximum(selectedNodeIdxes))
    newIdx = zero(UInt32)
    for oldIdx in selectedNodeIdxes
        newIdx += one(UInt32)
        oldIdx2newIdx[ oldIdx ] = newIdx
    end
    
    # since the root node have parent index 0, it can not be indexed directly
    @inbounds for (i, oldParentIdx) in enumerate(parents)
        if oldParentIdx > zero(UInt32)
            parents[i] = oldIdx2newIdx[oldParentIdx]
        end
    end
end