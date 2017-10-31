module NodeNets
include("DBFs.jl"); using .DBFs;
include("PointArrays.jl"); using .PointArrays;
using ..RealNeuralNetworks.SWCs

import LightGraphs
using Base.Cartesian
import BigArrays

const ZERO_UINT32 = convert(UInt32, 0)
const ONE_UINT32  = convert(UInt32, 1)
const OFFSET = (ZERO_UINT32, ZERO_UINT32, ZERO_UINT32)
const EXPANSION = (ONE_UINT32, ONE_UINT32, ONE_UINT32)

export NodeNet 

type NodeNet 
    # x,y,z,r
    nodeList               :: Vector{NTuple{4,Float32}}
    # connectivity matrix to represent edges
    # conn[2,3]=true means node 2 and 3 connected with each other
    # conn[3,2] should also be true since this is undirected graph
    connectivityMatrix  :: SparseMatrixCSC{Bool,UInt32}
end 

"""
    NodeNet( seg, obj_id; penalty_fn=alexs_penalty)
Perform the teasar algorithm on the passed binary array.
"""
function NodeNet{T}( seg::Array{T,3}; 
                        obj_id::T = convert(T,1), 
                        expansion::NTuple{3, UInt32} = EXPANSION,
                        penalty_fn::Function = alexs_penalty )
    # note that the object voxels are false and non-object voxels are true!
    # bin_im = DBFs.create_binary_image( seg, obj_id ) 
    points = PointArrays.from_seg(seg; obj_id=obj_id)
    NodeNet(points; expansion=expansion, penalty_fn=penalty_fn) 
end 

"""
    NodeNet(bin_im)
Parameters:
    bin_im: binary mask. the object voxel should be false, non-object voxel should be true
Return:
    nodeNet object
"""
function NodeNet(bin_im::Array{Bool,3}; offset::NTuple{3, UInt32} = OFFSET,
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
    NodeNet(points; dbf=dbf, penalty_fn=penalty_fn, expansion = expansion)
end 

"""
    NodeNet( points; penalty_fn = alexs_penalty )

  Perform the teasar algorithm on the passed Nxd array of points
"""
function NodeNet{T}( points::Array{T,2}; dbf=DBFs.compute_DBF(points),
                            penalty_fn::Function = alexs_penalty,
                            expansion::NTuple{3, UInt32} = EXPANSION)

    println("total number of points: $(size(points,1))")
  points, bbox_offset = shift_points_to_bbox( points );
  ind2node, max_dims = create_node_lookup( points );
  max_dims_arr = [max_dims...];#use this for rm_nodes, but ideally wouldn't
  sub2node = x -> ind2node[ sub2ind(max_dims, x[1],x[2],x[3]) ];#currently only used in line 48
    
  println("making graph (2 parts)");
  @time G, weights = make_neighbor_graph( points, ind2node, max_dims;)
  println("build dbf weights from penalty function ...")
  @time dbf_weights = penalty_fn( weights, dbf, G )

  #init
  #nonzeros SHOULD remove duplicates, but it doesn't so
  # I have to do something a bit more complicated
  _,nonzero_vals = findnz(ind2node);
  disconnected_nodes = IntSet( nonzero_vals );
  paths = Vector(); #holds array of nodeNet paths
  root_nodes = Vector{Int}(); #holds root nodes for each path
  destinations = Vector{Int}(); #host dest node for each path
  #set of nodes for which we've "inspected" already
  # removing their neighbors based on DBF
  inspected_nodes = Set{Int}();

  println("Finding paths")
  while length(disconnected_nodes) > 0

    root_ind = find_new_root_node( points[[disconnected_nodes...],:] );
    root_node = sub2node( root_ind );
    push!(root_nodes, root_node);
    println("root");
    println(root_node);

    #do a graph traversal to find farthest node and
    # find reachable nodes
    dsp_euclidean = LightGraphs.dijkstra_shortest_paths(G, root_node, weights);
    #and another to find the min DBF-weighted paths to all nodes
    dsp_dbf       = LightGraphs.dijkstra_shortest_paths(G, root_node, dbf_weights);

    #remove reachable nodes from the disconnected ones
    reachable_nodes = find(!(isinf(dsp_euclidean.dists)));
    #another necessary precaution for duplicate nodes
    reachable_nodes = intersect(reachable_nodes, disconnected_nodes);
    setdiff!(disconnected_nodes, reachable_nodes);
    empty!(inspected_nodes);

    while length(reachable_nodes) > 0

      #find the node farthest away from the root
      # by euc distance
      _, farthest_index = findmax( dsp_euclidean.dists[[reachable_nodes...]] );
      farthest_node = reachable_nodes[farthest_index];
      push!(destinations, farthest_node);
      println("dest")
      println(farthest_node)

      if farthest_node == root_node break end #this can happen apparently

      new_path = LightGraphs.enumerate_paths( dsp_dbf,farthest_node );

      push!(paths, new_path)
      #can't do this in-place with arrays
      #this fn call is getting ridiculous
      @time reachable_nodes = remove_path_from_rns( reachable_nodes, new_path,
                                              points, ind2node,
                                              dbf, max_dims_arr,
                                              inspected_nodes );

    end #while reachable nodes from root
  end #while disconnected nodes

  println("Consolidating Paths")
  path_nodes, path_edges = consolidate_paths( paths );
  node_radii = dbf[path_nodes];
  
    # build a new graph containing only the nodeNet nodes and edges
    nodes, edges = distill!(points, path_nodes, path_edges)

    conn = get_connectivity_matrix(edges)
    nodeList = Vector{NTuple{4,Float32}}()
    for i in 1:length(node_radii)
        push!(nodeList, (map(Float32,nodes[i,:])..., node_radii[i]))
    end 
    nodeNet = NodeNet(nodeList, conn)
    # add the offset from shift bounding box function
    @show bbox_offset
    add_offset!(nodeNet, bbox_offset)
    return nodeNet
end

function NodeNet(swc::SWC)
    nodeList = Vector{NTuple{4, Float32}}()
    connectivityMatrix = spzeros(Bool, length(swc), length(swc))
    for (index, point) in enumerate(swc) 
        push!(nodeList, (point.x, point.y, point.z, point.radius))
        if point.parent != -1
            # the connectivity matrix is symetric
            connectivityMatrix[index, point.parent] = true
            connectivityMatrix[point.parent, index] = true
        end 
    end 
    NodeNet(nodeList, connectivityMatrix)
end 

function SWCs.SWC(nodeNet::NodeNet)
    edges = NodeNets.get_edges(nodeNet)
    swc = SWC()
    sizehint!(swc, NodeNets.get_node_num(nodeNet))

    for node in NodeNets.get_node_list(nodeNet)
        point = SWCs.PointObj(0, node[1], node[2], node[3], node[4], -1)
        push!(swc, point)
    end
    # assign parents according to edge 
    for e in edges 
        swc[e[2]].parent = e[1]
    end  
    swc
end


##################### properties ###############################
function get_node_list(self::NodeNet) self.nodeList end 
function get_connectivity_matrix(self::NodeNet) self.connectivityMatrix end
function get_xyz(self::NodeNet) map(x->x[1:3], self.nodeList) end
function get_radii(self::NodeNet) map(x->x[4],  self.nodeList) end 
function get_node_num(self::NodeNet) length(self.nodeList) end
# the connectivity matrix is symmetric, so the connection is undirected
function get_edge_num(self::NodeNet) div(nnz(self.connectivityMatrix), 2) end

function get_edges(self::NodeNet) 
    edges = Vector{Tuple{UInt32,UInt32}}()
    conn = get_connectivity_matrix(self)
    I,J,V = findnz(conn)
    # only record the triangular part of the connectivity matrix
    for index in 1:length(I)
        if I[index] > J[index]
            push!(edges, (I[index], J[index]) )
        end 
    end 
    edges
end 

"""
    get_branch_point_num(self::NodeNet)

get number of branching points 
"""
function get_num_branch_point(self::NodeNet)
    conn = get_connectivity_matrix(self)
    num_branch_point = 0
    for col in size(conn, 2)
        if nnz(conn[:,col]) > 2
            num_branch_point += 1
        end
    end 
    return num_branch_point 
end

"""
assume that the graph is acyclic, no loop.
"""
function get_num_branches(self::NodeNet)
    get_num_branch_point(self) + 2
end

"""
    get_sholl_number(self::NodeNet, radius::AbstractFloat)
get the number of points which is in neurite and incounters with a sphere centered on root node 
"""
function get_sholl_number(self::NodeNet, radius::AbstractFloat; rootNodeIndex::Integer=1)
    shollNum = 0
    nodeList = get_node_list(self)
    rootNode = nodeList[ rootNodeIndex ]
    for edge in get_edges(self)
        node1 = nodeList[ edge[1] ]
        node2 = nodeList[ edge[2] ]
        r1 = distance_from( rootNode, node1 )
        r2 = distance_from( rootNode, node2 )
        if r1>radius != r2>radius
            shollNum += 1
        end 
    end 
    shollNum
end 

function Base.UnitRange(self::NodeNet)
    minCoordinates = [typemax(UInt32), typemax(UInt32), typemax(UInt32)]
    maxCoordinates = [ZERO_UINT32, ZERO_UINT32, ZERO_UINT32]
    for node in get_node_list(self)
        minCoordinates = map(min, minCoordinates, node[1:3])
        maxCoordinates = map(max, maxCoordinates, node[1:3])
    end 
    return [minCoordinates[1]:maxCoordinates[1], 
            minCoordinates[2]:maxCoordinates[2], 
            minCoordinates[3]:maxCoordinates[3]]
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
"""
function get_neuroglancer_precomputed(self::NodeNet)
    @show get_node_num(self)
    @show get_edge_num(self)
    # total number of bytes
    num_bytes = 4 + 4 + 4*3*get_node_num(self) + 4*2*get_edge_num(self)
    buffer = IOBuffer( num_bytes )
    # write the number of vertex, and edges
    write(buffer, UInt32(get_node_num(self)))
    write(buffer, UInt32(get_edge_num(self)))
    # write the node coordinates
    for node in get_node_list(self)
        write(buffer, [node[1:3]...])
    end
    # write the edges
    for edge in get_edges( self )
        # neuroglancer index is 0-based
        write(buffer, UInt32( edge[1]-ONE_UINT32 ))
        write(buffer, UInt32( edge[2]-ONE_UINT32 ))
    end
    bin = Vector{UInt8}(take!(buffer))
    close(buffer)
    return bin 
end 

"""
    get_connectivity_matrix(edges::Vector)

construct sparse connectivity matrix accordint to the edges 
"""
function get_connectivity_matrix(edges::Vector)
    I = [e[1] for e in edges]
    J = [e[2] for e in edges]
    return sparse([I...,J...],[J...,I...],ones(Bool, 2*length(edges)))
end 

###################### IO #################################
"""
    save(self::NodeNet, cellId::UInt32, d_json::Associative, d_bin::Associative)

save nodeNet in google cloud storage for neuroglancer visualization
the format is the same with meshes
"""
function save(self::NodeNet, cellId::Integer, d_bin::Associative)
    # get the bounding box of nodeNet and transfer to string representation
    # example string: 1432-1944_1264-1776_16400-16912
    rangeString = BigArrays.Indexes.unit_range2string( UnitRange(self) )
    # write the binary representation of nodeNet
    d_bin["$(cellId)"] = get_neuroglancer_precomputed(self)
end

"""
save binary file of point pairs
used in neuroglancer python interface to visualize the nodeNet 
"""
function save_edges(self::NodeNet, fileName::String)
    open(fileName, "w") do f
        nodeList = get_node_list(self)
        edges = get_edges(self)
        for edge in edges
            write(f, nodeList[ edge[1] ])
            write(f, nodeList[ edge[2] ])
        end 
    end
end 

##################### manipulate ############################
function add_offset!(self::NodeNet, offset::Union{Vector,Tuple} )
    @assert length(offset) == 3
    for i in 1:get_node_num(self)
        xyz = map((x,y)->x+y, self.nodeList[i][1:3],offset)
        self.nodeList[i] = (xyz..., self.nodeList[i][4])
    end 
end

function stretch_coordinates!(self::NodeNet, mip::Real)
    expansion = (2^mip, 2^mip, 1) 
    stretch_coordinates!(self, expansion)
end 
function stretch_coordinates!(self::NodeNet, expansion::Union{Vector, Tuple})
    @assert length(expansion) == 3
    nodeList = get_node_list(self)
    for (nodeIndex, node) in enumerate(nodeList)
        self.nodeList[ nodeIndex ] = (map((x,y)->x*y, node[1:3], expansion)..., 
                                      node[4]* (prod(expansion)^(1/3)))
    end 
end 

#---------------------------------------------------------------


"""

    shift_points_to_bbox( points )

  Normalize the point dimensions by subtracting the min
  across each dimension. This step isn't extremely necessary,
  but might be useful for compatibility with the MATLAB code.
  record the offset and add it back after building the nodeNet
"""
function shift_points_to_bbox( points )
  offset = minimum( points, 1 ) -1 ;
  # transform to 1d vector
  offset = vec( offset )
  @assert ndims(offset) == 1
  @assert length(offset) == 3
  @assert offset[3] < 20000
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
function create_node_lookup{T}( points::Array{T,2} )

  #need the max in order to define the bounds of that volume
  # tuple allows the point to be passed to fns as dimensions
  max_dims = ( maximum( points, 1 )... );
  max_dims = map(Int, max_dims)
  num_points = size(points,1);
  #creating the sparse vector mapping indices -> node ids
  ind2node = sparsevec( [sub2ind(max_dims, Vector{Int}(points[i,:])... ) for i=1:num_points ],
  	                    1:num_points,
  	                    prod(max_dims),
  	                    (x,y) -> x #ignore duplicates
  	                    )

  ind2node, max_dims
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
function make_neighbor_graph{T}( points::Array{T,2}, ind2node=nothing, max_dims=nothing;
                                expansion::NTuple{3, UInt32} = EXPANSION )

  if ind2node == nothing ind2node, max_dims = create_node_lookup(points) end

  #init
  num_points = size(points,1);
  g = LightGraphs.Graph(num_points);

  #26-connectivity neighborhood
  # weights computed by euc_dist to center voxel
  nhood = [[i,j,k] for i=1:3,j=1:3,k=1:3];
  map!( x-> (x .- [2,2,2]).*[expansion ...] , nhood );
  nhood_weights = map( norm, nhood );

  #only adding weights for non-duplicate nodes
  _,non_duplicates = findnz( ind2node );

  #init - descriptors of the edges for graph
  # weights later
  is = Vector{Int}();
  js = Vector{Int}();
  ws = Vector{Float64}();

  dest_sub = Array{Int}((3,))

  for i in eachindex(nhood)

  	#only need to explicitly tell LightGraphs to add the edges in
  	# one direction (undirected), though doing this properly
    # can be tricky
    direction = ind2sub((3,3,3),i)
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
      dest_ind = sub2ind( max_dims, dest_sub[1],dest_sub[2],dest_sub[3] );
      dest_node = ind2node[dest_ind];

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
  M = maximum(dbf,1)*1.01

  p_v = (1 - dbf./M ).^(16);

  #DBF penalty is defined by the destination node, only need js
  is,js = findn(weights);
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

    find_new_root_node( points )

  Extracts the point with the lowest linear index from the
  Array of passed points
"""
function find_new_root_node( points )

  res = points;
  for i in 3:-1:1
    res = filter_rows_by_min( res, i );
  end

  res
end


function filter_rows_by_min{T}( arr::Array{T,2}, dim )

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

    remove_path_from_rns( reachable_nodes, path, points, sub2node, dbf, max_dims, scale_param, const_param)

  Identifies the nodes to remove from the reachable nodes within the graph.
  Probably the ugliest function here, and should be optimized later

  TO OPTIMIZE
"""
function remove_path_from_rns( reachable_nodes::Vector, path::Vector,
  points, ind2node, dbf, max_dims, inspected_nodes,
  scale_param::Int=2, const_param::Int=4 );

  r = dbf[path]*scale_param + const_param;

  to_remove = Set{Int}();
  for i in eachindex(path)
    node = path[i];
    if !(node in inspected_nodes)
        push!(inspected_nodes, node);
	    union!( to_remove, nodes_within_radius( points[node,:], ind2node, r[i], max_dims ) );
    end
  end

  return setdiff(reachable_nodes, to_remove);
end


"""

    nodes_within_radius( sub, sub2node, r, max_dims )

  Identifies the node indices within the subscript
"""
function nodes_within_radius{T}( sub::Array{T,1}, ind2node, r, max_dims::Vector );

  beginning = convert(Vector{Int}, ceil(max(sub[:] .- r,1)));
  ending    = convert(Vector{Int}, floor(min(sub[:] .+ r, max_dims)));
  ind::Int = convert(Int,0);
    max_dims = Vector{Int}(max_dims)

  nodes = Set{T}();
  for x in beginning[1]:ending[1]
    for y in beginning[2]:ending[2]
      for z in beginning[3]:ending[3]
        ind = x + (y-1)*max_dims[1] + (z-1)*max_dims[2]*max_dims[1];
        push!( nodes, ind2node[ind] );
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

    nodes = Set{Int}();
    edges = Set{NTuple{2,Int}}();

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
        @assert path_length > 0

        for i in 1:(path_length-1)
            push!( edges, (path[i],path[i+1]) );
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
function distill!{T}(point_array::Array{T,2}, 
                     path_nodes::Vector{Int}, path_edges::Vector) 
    num_nodes = length(path_nodes) 
    num_edges = length(path_edges)
    nodes = Array(T, (num_nodes, 3))
    # map the old path node id to new id 
    id_map = Dict{T, T}()
    sizehint!(id_map, num_nodes)
    # build new nodes and the id map
    for i in 1:num_nodes 
        nodes[i,:] = point_array[path_nodes[i], :]
        id_map[path_nodes[i]] = i
    end 

    # rebuild the edges
    for i in 1:num_edges
        path_edges[i] = (id_map[path_edges[i][1]], id_map[path_edges[i][2]])
    end
    return nodes, path_edges
end 

end # module
