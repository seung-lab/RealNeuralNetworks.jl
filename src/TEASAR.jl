export teasar 

import LightGraphs
import DataStructures: IntSet 

include("DBFs.jl"); using .DBFs;

# control the removing points around path based on DBF
const REMOVE_PATH_SCALE = 3 
const REMOVE_PATH_CONST = 4

"""
    teasar( seg, obj_id; penalty_fn=alexs_penalty)
Perform the teasar algorithm on the passed binary array.
"""
function teasar( seg::Array{T,3}; 
                     obj_id::T = convert(T,1), 
                     expansion::NTuple{3, UInt32} = EXPANSION,
                     penalty_fn::Function = alexs_penalty ) where T
    # note that the object voxels are false and non-object voxels are true!
    # bin_im = DBFs.create_binary_image( seg, obj_id ) 
    points = from_seg(seg; obj_id=obj_id)
    teasar(points; expansion=expansion, penalty_fn=penalty_fn) 
end 

"""
    teasar(bin_im)
Parameters:
    bin_im: binary mask. the object voxel should be false, non-object voxel should be true
Return:
    nodeNet object
"""
function teasar(bin_im::Union{BitArray, Array{Bool,3}}; 
                 offset::NTuple{3, UInt32} = OFFSET,
                 expansion::NTuple{3, UInt32} = EXPANSION,
                 penalty_fn::Function = alexs_penalty)
        # transform segmentation to points
    points = from_binary_image(bin_im)
    
    println("computing DBF");
    # boundary_point_indexes = PointArrays.get_boundary_point_indexes(points, seg; obj_id=obj_id)
    #@time DBF = DBFs.compute_DBF( points, boundary_point_indexes );
    @time DBF = DBFs.compute_DBF(points)
    # @time dbf = DBFs.compute_DBF(points, bin_im)

    add_offset!(points, offset)
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
    nodeNet = distill!(points, path_nodes, path_edges)
    set_radii!(nodeNet, node_radii)

    # add the offset from shift bounding box function
    bbox_offset = map(Float32, bbox_offset)
    add_offset!(nodeNet, bbox_offset)
    nodeNet
end


"""

    translate_to_origin!( points::Matrix{T} )

  Normalize the point dimensions by subtracting the min
  across each dimension. This step isn't extremely necessary,
  but might be useful for compatibility with the MATLAB code.
  record the offset and add it back after building the nodeNet
"""
function translate_to_origin!( points::Matrix{T} ) where T
    offset = minimum( points, dims=1 ) .- one(T) ;
    # transform to 1d vector
    offset = vec( offset )
    @assert length(offset) == 3
    @assert offset[3] < T(20000)
    points[:,1] .-= offset[1]
    points[:,2] .-= offset[2]
    points[:,3] .-= offset[3]
  
    points, offset
end 

"""

    create_node_lookup( points )

  Abstractly represent the points as a volume by means of
  linear indexes into a sparse vector.
"""
function create_node_lookup( points::Matrix{T} ) where T

  #need the max in order to define the bounds of that volume
  # tuple allows the point to be passed to fns as dimensions
  max_dims = ( maximum( points, dims=1 )... ,);
  max_dims = map(Int, max_dims)
  num_points = size(points,1);
  #creating the sparse vector mapping indices -> node ids
  volumeIndex2NodeIndex = sparsevec( [
                        (LinearIndices(max_dims))[Vector{UInt32}(points[i,:])...,] 
                            for i=1:num_points ], 1:num_points, prod(max_dims), (x,y) -> x )

  volumeIndex2NodeIndex, max_dims
end

"""

    make_neighbor_graph( points )

  Creates a LightGraphs graph from a point cloud, such that
  each point has a node (indexed in the same order as the rows),
  and each node is connected to its neighbors (by 26-connectivity).

  Also returns a sparse matrix of weights for these edges equal to the
  euclidean distance between the indices of each point. These can
  be weighted or modified easily upon return.
"""
function make_neighbor_graph( points::Array{T,2}, volumeIndex2NodeIndex=nothing, max_dims=nothing;
                             expansion::NTuple{3, UInt32} = EXPANSION ) where T

  if volumeIndex2NodeIndex == nothing volumeIndex2NodeIndex, max_dims = create_node_lookup(points) end

  #init
  num_points = size(points,1);
  g = LightGraphs.Graph(num_points);

  #26-connectivity neighborhood
  # weights computed by euc_dist to center voxel
  nhood = [[i,j,k] for i=1:3,j=1:3,k=1:3];

  map!( x-> (x .- [2,2,2]).*[expansion ...] , nhood, nhood)
  nhood_weights = map( norm, nhood );

  #only adding weights for non-duplicate nodes
  _,non_duplicates = findnz( volumeIndex2NodeIndex );

  #init - descriptors of the edges for graph
  # weights later
  is = Vector{Int}();
  js = Vector{Int}();
  ws = Vector{Float64}();

  dest_sub = Array{Int}(undef, 3)

  for i in eachindex(nhood)

    #only need to explicitly tell LightGraphs to add the edges in
    # one direction (undirected), though doing this properly
    # can be tricky
    #direction = ind2sub((3,3,3),i)
    direction = Tuple(CartesianIndices((3,3,3))[i]) 

    if (direction[1] == 3 ||
       (direction[1] == 2 && 2*direction[2] + direction[3] >= 6)) continue end

    for p in non_duplicates

      #location of potential node to connect
      for j=1:3 dest_sub[j] = nhood[i][j] + points[p,j]' end

      #bounds_checking
      if (dest_sub[1] < 1 ||
          dest_sub[2] < 1 ||
          dest_sub[3] < 1 ) continue end

      if (dest_sub[1] > max_dims[1] ||
          dest_sub[2] > max_dims[2] ||
          dest_sub[3] > max_dims[3] ) continue end

      #fetching the id of the proposed node location
      dest_ind = (LinearIndices(max_dims))[dest_sub[1:3]...,]
      dest_node = volumeIndex2NodeIndex[dest_ind];

      if dest_node != 0 #valid edge
        LightGraphs.add_edge!(g,p,dest_node)
        push!(is,p); push!(js,dest_node); push!(ws,nhood_weights[i]);
      end

    end #for p
  end #for i in nhood

  #creating weight matrix
  g_weights = sparse(vcat(is,js),vcat(js,is),vcat(ws,ws), num_points, num_points);

  g,g_weights
end


"""

    alexs_penalty( weights, dbf )

  Returns a version of the edge weights, modified by the DBF-based
  teasar penalty:

  w = w * 5000 .* (1 - dbf/maxDBF).^(16)

  The factor of 5000 shouldn't make a difference, but I've
  left it here
"""
function alexs_penalty( weights, dbf, G, mult=true )

  # *1.01 ensures later quotient is >0
  M = maximum(dbf, dims=1)*1.01

  p_v = (1 .- dbf./M ) .^ 16;

  #DBF penalty is defined by the destination node, only need js
  is, js = begin 
      I = findall(!iszero, weights) 
      (getindex.(I, 1), getindex.(I, 2))
  end 
  ws = nonzeros(weights);

  if mult
    new_ws = ws .* p_v[js];
  else
    new_ws = ws .+ p_v[js];
  end

  sparse( is, js, new_ws, size(weights)... );
end


"""

    local_max_multiplicative_penalty( weights, dbf, G )

  Returns a version of the edge weights, modified by the DBF-based
  penalty

  w = w * (1 - DBFdest/DBFstar)

  Where DBFstar is the maximum DBF value of an outwards neighbor for
  each node
"""
function local_max_multiplicative_penalty( weights, dbf, G )
    ns = LightGraphs.vertices(G)
    dbf_star = Vector{Float64}( length(ns) );

    #this is ok since ns is always an int range 1:N
    for n in ns
        neighborhood = LightGraphs.out_neighbors(G, n);

        if length(neighborhood) > 0
            dbf_star[n] = maximum(dbf[neighborhood]);
        else
            dbf_star[n] = 1
        end
    end

    is, js = findn(weights)
    ws = nonzeros(weights)
    new_ws = Vector{Float64}(length(ws));

    for i in eachindex(is)
        dbf_dest = dbf[js[i]];
        dbf_s = dbf_star[is[i]];

        new_ws[i] = ws[i] * (1 - dbf_dest/dbf_s);
    end

    sparse( is, js, new_ws, size(weights)... );
end

"""

    find_new_root_node( dbf )

  Extracts the point with the largest dbf from the 
  Array of passed points 
"""
@inline function find_new_root_node_id(dbf::DBF, disconnectedNodeIdSet::IntSet)
    disconnectedNodeIndexList = [disconnectedNodeIdSet...]
    disconnectedNodeDBFList = dbf[disconnectedNodeIndexList]
    rootIndexInDisconnectedList = argmax( disconnectedNodeDBFList )
    rootIndex = disconnectedNodeIndexList[rootIndexInDisconnectedList]
    return rootIndex
end 

"""

    find_new_root_node( points )

  Extracts the point with the lowest linear index from the
  Array of passed points 
"""
@inline function find_new_root_node_V1( points::Array{T,2} ) where T

    res = points;
    for i in 3:-1:1
        res = filter_rows_by_min( res, i );
    end

    res
end

function filter_rows_by_min( arr::Array{T,2}, dim ) where T
    selected_rows = arr[:,dim] .== minimum(arr,1)[dim];
    arr[selected_rows,:];
end 

function create_dummy_matrix( s )
    points = hcat( findn(ones(s,s,s))... );
    i = div(size(points,1),2)+1;

    points = points[[1:i-1;i+1:end],:];
    points;
end

"""

    remove_path_from_rns!( reachableNodeList, path, points, sub2node, dbf, max_dims, scale_param, const_param)

  Identifies the nodes to remove from the reachable nodes within the graph.
  Probably the ugliest function here, and should be optimized later

  TO OPTIMIZE
"""
function remove_path_from_rns!( reachableNodeList::Vector, path::Vector{Int},
                points, volumeIndex2NodeIndex, dbf, max_dims, inspectedNodeIdList,
                scale_param::Integer=REMOVE_PATH_SCALE, 
                const_param::Integer=REMOVE_PATH_CONST )
    radii = dbf[path].*scale_param .+ const_param;

    to_remove = Set{Int}();
    for i in eachindex(path)
        nodeIndex = path[i];
        if !(nodeIndex in inspectedNodeIdList)
            push!(inspectedNodeIdList, nodeIndex);
            union!( to_remove, nodes_within_radius(points[nodeIndex,:], 
                    volumeIndex2NodeIndex, radii[i], max_dims ) );
        end
    end

    return setdiff(reachableNodeList, to_remove);
end


"""

    nodes_within_radius( sub, sub2node, r, max_dims )

  Identifies the node indices within the subscript
"""
function nodes_within_radius( sub::Array{T,1}, volumeIndex2NodeIndex, r, max_dims::Vector ) where T;

    beginning = convert(Vector{Int}, ceil.(max.(sub[:] .- r,1)));
    ending    = convert(Vector{Int}, floor.(min.(sub[:] .+ r, max_dims)));
    ind::Int = convert(Int,0);
      max_dims = Vector{Int}(max_dims)

    nodes = Set{T}();
    for x in beginning[1]:ending[1]
        for y in beginning[2]:ending[2]
            for z in beginning[3]:ending[3]
                ind = x + (y-1)*max_dims[1] + (z-1)*max_dims[2]*max_dims[1];
                push!( nodes, volumeIndex2NodeIndex[ind] );
                #don't need to worry about the '0 node' here as
                # the only result that matters is the setdiff from
                # the reachable nodes
            end
        end
    end

    nodes
end


"""

    consolidate_paths( path_list )

  Extracts the unique nodes and edges among all paths in the list
"""
function consolidate_paths( path_list::Vector )

    nodes = Set{UInt32}();
    edges = Set{NTuple{2,UInt32}}();

    # estimate the total number of nodes and edges
    num_nodes = 0
    num_edges = 0
    for path in path_list
        path_length = length(path)
        num_nodes += path_length 
        num_edges += path_length-1 
    end 
    # give some hints for memory allocation
    sizehint!(nodes, num_nodes)
    sizehint!(edges, num_edges)

    for path in path_list
        path_length = length(path);
        @assert path_length > 1

        for i in 1:(path_length-1)
            @assert path[i] != path[i+1]
            push!( edges, (path[i], path[i+1]) );
            push!( nodes, path[i] );
        end
        push!( nodes, path[path_length] );
    end

    collect(nodes), collect(edges)
end

"""
Parameters:
===========
    point_array: a Nx3 array recording all the voxel coordinates inside the object.
    path_nodes: the indexes of nodeNet voxels in the point_array 
    path_edges: the index pairs of nodeNet in the point_array  

In the path_nodes, the nodes were encoded as the index of the point_array. Since we do 
not need point_array anymore, which could be pretty big, we'll only reserve the nodeNet 
coordinates and let the gc release the point_array. 
the same applys to path_edges 
"""
function distill!(point_array::Array{T,2}, 
                  path_nodes::Vector, path_edges::Vector) where T
    @assert size(point_array, 2) == 3
    nodeNum = length(path_nodes) 
    edgeNum = length(path_edges)
    nodes = Matrix{Float32}(undef, 4, nodeNum)

    # map the old path node id to new id
    # the old node index number is more than new nodes
    oldIdx2newIdx = zeros(UInt32, size(point_array, 1)) 
    # build new nodes and the id map
    for i in one(UInt32):UInt32(nodeNum)
        # we'll deal with radius later 
        nodes[1:3, i] = Float32.(point_array[path_nodes[i], :])
        oldIdx2newIdx[path_nodes[i]] = i
    end 

    # rebuild the edges
    childIdxes = zeros(UInt32, edgeNum)
    parentIdxes = zeros(UInt32, edgeNum)
    for i in 1:edgeNum
        childIdx, parentIdx = path_edges[i]
        # map to new index 
        childIdx = oldIdx2newIdx[ childIdx ]
        parentIdx = oldIdx2newIdx[ parentIdx ]
        @assert childIdx != parentIdx
        childIdxes[i] = childIdx
        parentIdxes[i] = parentIdx
    end

    connectivityMatrix = sparse(childIdxes, parentIdxes, true, nodeNum, nodeNum)

    NodeNet(nodes, connectivityMatrix)
end 