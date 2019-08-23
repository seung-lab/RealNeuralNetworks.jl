module NodeNets
include("TEASAR.jl")

import LinearAlgebra: norm
import DelimitedFiles: readdlm, writedlm
const ONE_UINT32 = UInt32(1)
using SparseArrays 
using NearestNeighbors

export NodeNet 


const OFFSET = (zero(UInt32), zero(UInt32), zero(UInt32))
# rescale the skeleton
const EXPANSION = (one(UInt32), one(UInt32), one(UInt32))

mutable struct NodeNet{T}  
    # x,y,z,r
    classes             :: Vector{UInt8}
    nodeArray           :: Matrix{T}
    # connectivity matrix to represent edges
    # conn[2,3]=true means node 2's parent is node 3 
    connectivityMatrix  :: SparseMatrixCSC{Bool,UInt32}
end 

function NodeNet( data::Matrix{T} ) where T
    @assert size(data, 2) == 7
    data = data[sortperm(data[:, 1]), :]
    classes = UInt8.(data[:, 2])
    nodeArray = Float32.(data[:, 3:6])
    parents = Int32.(data[:, 7])
    connectivityMatrix = parents_to_connectivity_matrix(parents)
    NodeNet(classes, nodeArray, connectivityMatrix)
end

@inline function Matrix{T}(self::NodeNet) where T
    parents = connectivity_matrix_to_parents(self.connectivityMatrix)
    hcat(self.classes, self.nodeArray, parents)
end

function NodeNet(nodeArray::Matrix{T}, 
                    connectivityMatrix::SparseMatrixCSC{Bool, UInt32}) where T
    @assert size(nodeArray, 2) == 4
    classes = zeros(UInt8, size(nodeArray, 1))
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
    teasar( points; penalty_fn = alexs_penalty )

  Perform the teasar algorithm on the passed Nxd array of points
"""
function teasar( points::Matrix{T}; dbf::DBF=DBFs.compute_DBF(points),
                         penalty_fn::Function = alexs_penalty,
                         expansion::NTuple{3, UInt32} = EXPANSION) where T
    @assert length(dbf) == size(points, 1)
    println("total number of points: $(size(points,1))")
    points, bbox_offset = translate_to_origin!( points );
    # volumeIndex2NodeIndex represent points as a sparse vector containing the whole volume
    # in this way, we can fetch the node directly according to the coordinate.
    # the coordinate should be transfered to vector index though
    volumeIndex2NodeId, max_dims = create_node_lookup( points );
    max_dims_arr = [max_dims...];#use this for rm_nodes, but ideally wouldn't
    # transfer the coordinate to node 
    # sub2node = x -> volumeIndex2NodeIndex[ sub2ind(max_dims, x[1],x[2],x[3]) ];#currently only used in line 48
    
    println("making graph (2 parts)");
    @time G, weights = make_neighbor_graph( points, volumeIndex2NodeId, max_dims;)
    println("build dbf weights from penalty function ...")
    @time dbf_weights = penalty_fn( weights, dbf, G )

    #init
    #nonzeros SHOULD remove duplicates, but it doesn't so
    # I have to do something a bit more complicated
    _,nonzero_vals = findnz(volumeIndex2NodeId);
    disconnectedNodeIdSet = IntSet( nonzero_vals );
    pathList = Vector(); # holds vector of nodeNet paths
    destinationNodeIdList = Vector{Int}(); #host dest node for each path
    # set of nodes for which we've "inspected" already
    # removing their neighbors based on DBF
    inspectedNodeIdList = Set{Int}();

    println("Finding paths")
    while length(disconnectedNodeIdSet) > 0
        rootNodeId = find_new_root_node_id( dbf, disconnectedNodeIdSet );
        @assert rootNodeId in disconnectedNodeIdSet 

        #do a graph traversal to find farthest node and
        # find reachable nodes
        dsp_euclidean = LightGraphs.dijkstra_shortest_paths(G, rootNodeId, weights);
        #and another to find the min DBF-weighted paths to all nodes
        dsp_dbf       = LightGraphs.dijkstra_shortest_paths(G, rootNodeId, dbf_weights);

        #remove reachable nodes from the disconnected ones
        reachableNodeIdList = findall(.!(isinf.(dsp_euclidean.dists)));
        #another necessary precaution for duplicate nodes
        reachableNodeIdList = intersect(reachableNodeIdList, disconnectedNodeIdSet);
        setdiff!(disconnectedNodeIdSet, reachableNodeIdList);
        empty!(inspectedNodeIdList);

        while length(reachableNodeIdList) > 0

            #find the node farthest away from the root
            # by euc distance
            _, farthestNodeIndex = findmax( dsp_euclidean.dists[[reachableNodeIdList...]] );
            farthestNodeId = reachableNodeIdList[farthestNodeIndex];
            push!(destinationNodeIdList, farthestNodeId);
            println("dest node index: $(farthestNodeId)")

            if farthestNodeId == rootNodeId break end #this can happen apparently

            new_path = LightGraphs.enumerate_paths( dsp_dbf, farthestNodeId );

            push!(pathList, new_path)
            #can't do this in-place with arrays
            #this fn call is getting ridiculous
            @time reachableNodeIdList = remove_path_from_rns!( reachableNodeIdList, 
                                                        new_path, points, volumeIndex2NodeId,
                                                        dbf, max_dims_arr,
                                                        inspectedNodeIdList );
        end #while reachable nodes from root
    end #while disconnected nodes

    println("Consolidating Paths...")
    path_nodes, path_edges = consolidate_paths( pathList );
    node_radii = dbf[path_nodes];

    # build a new graph containing only the nodeNet nodes and edges
    nodes, edges = distill!(points, path_nodes, path_edges)

    conn = get_connectivity_matrix(edges)
    nodeList = Vector{NTuple{4,Float32}}()
    sizehint!(nodeList, length(node_radii))
    for i in 1:length(node_radii)
        push!(nodeList, (map(Float32,nodes[i,:])..., node_radii[i]))
    end

    nodeArray = node_list_to_array(nodeList)
    nodeNet = NodeNet(nodeArray, conn)
    # add the offset from shift bounding box function
    bbox_offset = map(Float32, bbox_offset)
    @show bbox_offset
    add_offset!(nodeNet, bbox_offset)
    return nodeNet
end

function parents_to_connectivity_matrix(parents::Vector{Int32})
    nodeNum = length(parents)

    childNodeIdxList = Vector{UInt32}()
    sizehint!(childNodeIdxList, nodeNum)
    parentNodeIdxList = Vector{UInt32}()
    sizehint!(parentNodeIdxList, nodeNum)

    rootNodeParent = -one(Int32)

    @inbounds for childNodeIdx in 1:nodeNum
        parentNodeIdx = parents[childNodeIdx]
        if parentNodeIdx != rootNodeParent
            push!(childNodeIdxList, childNodeIdx)
            push!(parentNodeIdxList, parentNodeIdx)
        end
    end

    connectivityMatrix = sparse(childNodeIdxList, parentNodeIdxList, 
                                                true, nodeNum, nodeNum)
    connectivityMatrix
end

function connectivity_matrix_to_edges(conn::SparseMatrixCSC{Bool, UInt32})
    I,J,_ = findnz(conn)
    edgeNum = length(I)
    # the first one is child, and the second one is parent
    edges = Matrix{UInt32}(undef, edgeNum, 2)
    # only record the triangular part of the connectivity matrix
    for index in 1:edgeNum
        edges[index, 1] = I[index]
        edges[index, 2] = J[index]
    end 
    edges
end 

function connectivity_matrix_to_parents(connectivityMatrix::SparseMatrixCSC{Bool, UInt32})
    nodeNum = size(connectivityMatrix, 1)
    parents = ones(Int32, nodeNum) .* -one(Int32)
    
    edges = connectivity_matrix_to_edges(connectivityMatrix)
    @inbounds for i in 1:size(edges, 1)
        parents[edges[i, 1]] = edges[i, 2]
    end  
    parents
end

##################### properties ###############################
@inline function get_node_array(self::NodeNet) self.nodeArray end 

@inline function get_node_list(self::NodeNet{T}) where T 
    nodeArray = get_node_array(self)
    N = size(nodeArray, 1)
    nodeList = Vector{NTuple{4,T}}()
    sizehint!(nodeList, N)
    for i in 1:N 
        node = tuple(nodeArray[i, :]...)
        push!(nodeList, node)
    end
    nodeList
end

@inline function get_connectivity_matrix(self::NodeNet) self.connectivityMatrix end
@inline function get_classes(self::NodeNet) self.classes end 
@inline function get_radii(self::NodeNet) self.nodeArray[:, 4] end 
function get_node_num(self::NodeNet; class::Union{Nothing, UInt8}=nothing)
    if class == nothing 
        return length(self.classes) 
    else 
        return count(self.classes .== class)
    end 
end

@inline function get_parents(self)
    connectivityMatrix = get_connectivity_matrix(self)
    connectivity_matrix_to_parents(connectivityMatrix)
end

"""
    node_list_to_array(nodeList::Vector{NTuple{4,T}}) where T

transform a list of nodes to array. The size of array is (4, N).
The first axis is the x,y,z,r, the second axis is the nodes.
"""
function node_list_to_array(nodeList::Vector{NTuple{4,T}}) where T
    N = length(nodeList)
    ret = Array{T,2}(undef, 4, N)
    @inbounds for (i,node) in enumerate(nodeList)
        ret[:,i] = [node...]
    end
    ret
end

""" 
the connectivity matrix is symmetric, so the connection is undirected
"""
@inline function get_edge_num(self::NodeNet) div(nnz(self.connectivityMatrix), 2) end

@inline function get_edges(self::NodeNet) 
    conn = get_connectivity_matrix(self)
    connectivity_matrix_to_edges(conn)
end 

@inline function is_terminal_node(self::NodeNet, i::Integer)
    connectivityMatrix = get_connectivity_matrix(self)
    dropzeros!(connectivityMatrix)
    return nnz(connectivityMatrix[i,:]) < 2
end

function get_terminal_node_id_list(self::NodeNet)
    terminalNodeIdList = Vector{Int}()
    connectivityMatrix = get_connectivity_matrix(self)
    dropzeros!(connectivityMatrix)
    for i in 1:get_node_num(self)
        if nnz(connectivityMatrix[i,:]) < 2
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
        if nnz(connectivityMatrix[i,:]) > 2
            push!(branchingNodeIdList, i)
        end 
    end
    return branchingNodeIdList
end 

"""
    get_segment_point_num(self::NodeNet)

get number of branching points 
"""
@inline function get_branching_node_num(self::NodeNet)
    get_branching_node_id_list(self) |> length
end

"""
assume that the graph is acyclic, no loop.
"""
@inline function get_num_segmentes(self::NodeNet)
    get_branching_node_num(self) + 2
end

@inline function euclidean_distance(node1::Tuple, node2::Tuple)
    norm([node1[1:3]...] .- [node2[1:3]...])
end

"""
    get_sholl_number(self::NodeNet, radius::AbstractFloat)
get the number of points which is in neurite and incounters with a sphere centered on root node 
"""
function get_sholl_number(self::NodeNet, radius::AbstractFloat; rootNodeIndex::Integer=1)
    error("depracated, use the function in Neurons.")
    shollNum = 0
    nodeList = get_node_list(self)
    rootNode = nodeList[ rootNodeIndex ]
    for edge in get_edges(self)
        node1 = nodeList[ edge[1] ]
        node2 = nodeList[ edge[2] ]
        r1 = euclidean_distance( rootNode, node1 )
        r2 = euclidean_distance( rootNode, node2 )
        if r1>radius != r2>radius
            shollNum += 1
        end 
    end 
    shollNum
end 

#################### Setters ############################################
function set_radius!(self::NodeNet, radius::Float32)
    nodeArray = get_node_array(self)
    nodeArray[:, 4] .= radius
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

@inline function Base.getindex(self::NodeNet, i::Integer)
    self.nodeArray[i, :]
end 

function Base.UnitRange(self::NodeNet)
    minCoordinates = [typemax(UInt32), typemax(UInt32), typemax(UInt32)]
    maxCoordinates = [zero(UInt32), zero(UInt32), zero(UInt32)]
    for i in 1:length(self)
        node = self[i]
        minCoordinates = map(min, minCoordinates, node[1:3])
        maxCoordinates = map(max, maxCoordinates, node[1:3])
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
    @assert N>=3 && N<=4
    find_closest_node_id(self, [point[1:3]...])
end

function find_closest_node_id(self::NodeNet{T}, point::Vector{T}) where T
    nodeArray = get_node_array(self)
    kdtree = KDTree(nodeArray' |> Matrix; leafsize=10)
    idxs, _ = knn(kdtree, point, 1)
    return idxs[1]
end

"""
    get_total_path_length( self::NodeNet )
accumulate all the euclidean distance of edges 
"""
function get_total_path_length( self::NodeNet{T} ) where T
    nodeArray = get_node_array(self)
    edges = get_edges(self)
    
    totalPathLength = zero(T)
    @inbounds for edgeIdx in 1:size(edges, 2)
        nodeIdx1 = edges[edgeIdx, 1]
        nodeIdx2 = edges[edgeIdx, 2]
        node1 = view(nodeArray, nodeIdx1, 1:3)
        node2 = view(nodeArray, nodeIdx2, 1:3)
        totalPathLength += norm(node1 .- node2)
    end
    totalPathLength
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
    @show get_node_num(self)
    @show get_edge_num(self)
    # total number of bytes
    num_bytes = 4 + 4 + 4*3*get_node_num(self) + 4*2*get_edge_num(self)
    io = IOBuffer( read=false, write=true, maxsize=num_bytes )
    
    # write the number of vertex, and edges
    write(io, UInt32(get_node_num(self)))
    write(io, UInt32(get_edge_num(self)))
    
    # write the node coordinates
    nodeArray = get_node_array(self)
    write(io, nodeArray[1:3, :]' |> Matrix)
    
    # write the edges
    edges = get_edges(self)
    # neuroglancer index start from 0
    edges = UInt32.(edges) .- one(UInt32)
    write(io, edges' |> Matrix)
    
    bin = Vector{UInt8}(take!(io))
    close(io)
    return bin 
end 

function deserialize(data::Vector{UInt8})
    # a pointObj is 21 byte
    @assert mod(length(data), 21) == 0 "the binary file do not match the byte layout of pointObj."
    nodeNum = div(length(data), 21)
    classes = Vector{UInt8}(undef, nodeNum)
    nodeArray = Matrix{Float32}(undef, nodeNum, 4)
    parents = Vector{Int32}(undef, nodeNum) 

    @inbounds for i in 1:nodeNum
        nodeData = view(data, (i-1)*21+1 : i*21)
        classes[i] = nodeData[1]
        nodeArray[i, :] = reinterpret(Float32, nodeData[2:17])
        parents[i] = reinterpret(Int32, nodeData[18:21])[1] 
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
        write(io, nodeArray[i, :])
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
    NodeNet( data )    
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
    xs = view(truncedNodeArray, :, 1)
    ys = view(truncedNodeArray, :, 2)
    zs = view(truncedNodeArray, :, 3)
    rs = view(truncedNodeArray, :, 4) 

    @show size(nodeArray), size(truncedNodeArray)

    nodeNum = length(classes)
    data = zip(1:nodeNum, classes, xs, ys, zs, rs, parents)
    writedlm(file_name, data, ' ', dims=(nodeNum, 7))
end

##################### manipulate ############################
@inline function add_offset!(self::NodeNet{T}, offset::NTuple{3,T} ) where T
    add_offset!(self, [offset...])
end

function add_offset!(self::NodeNet{T}, offset::Vector{T} ) where T
    @assert length(offset) == 3
    nodeArray = get_node_array(self)
    @inbounds for i in size(nodeArray, 1)
        nodeArray[i, 1:3] .+= offset
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
    @inbounds for i in 1:size(nodeArray, 1)
        nodeArray[i, :] .*= expansion
    end
end 

#---------------------------------------------------------------
end # module
