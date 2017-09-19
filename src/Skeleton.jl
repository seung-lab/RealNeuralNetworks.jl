module Skeleton

import LightGraphs
import ..BWDists 
import ..Points 
using Base.Cartesian

export skeletonize 

#---------------------------------------------------------------
"""
    skeletonize( seg; penalty_fn=alexs_penalty)
Perform the teasar algorithm on the passed binary array.
"""
function skeletonize{T}( seg::Array{T,3}; obj_id::T=convert(T,1), 
                   penalty_fn=alexs_penalty )
   # transform segmentation to points
   points = Points.from_seg(seg; obj_id=obj_id)
   skeletonize(points; penalty_fn=penalty_fn)
end 

"""
    skeletonize( points, scale_param, const_param, source_points )

  Perform the teasar algorithm on the passed Nxd array of points
"""
function skeletonize{T}( points::Array{T,2}; 
                   penalty_fn=alexs_penalty )#, scale_param, const_param, source_points )

  points = shift_points_to_bbox( points );
  ind2node, max_dims = create_node_lookup( points );
  max_dims_arr = [max_dims...];#use this for rm_nodes, but ideally wouldn't
  sub2node = x -> ind2node[ sub2ind(max_dims, x[1],x[2],x[3]) ];#currently only used in line 48

  println("computing DBF");
  @time DBF = compute_DBF( points );

  println("making graph (2 parts)");
  @time G, weights = make_neighbor_graph( points, ind2node, max_dims );
  @time dbf_weights = penalty_fn( weights, DBF, G );

  #init
  #nonzeros SHOULD remove duplicates, but it doesn't so
  # I have to do something a bit more complicated
  _,_,nonzero_vals = findnz(ind2node);
  disconnected_nodes = IntSet( nonzero_vals );
  paths = Vector(); #holds array of skeleton paths
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
                                              DBF, max_dims_arr,
                                              inspected_nodes );

    end #while reachable nodes from root
  end #while disconnected nodes

  println("Consolidating Paths")
  path_nodes, path_edges = consolidate_paths( paths );
  node_radii = DBF[path_nodes];

  path_edges, path_nodes, root_nodes, node_radii, destinations
end
#---------------------------------------------------------------


"""

    shift_points_to_bbox( points )

  Normalize the point dimensions by subtracting the min
  across each dimension. This step isn't extremely necessary,
  but might be useful for compatibility with the MATLAB code.
"""
function shift_points_to_bbox( points )
  min_dims = minimum( points, 1 );
  points .-= (min_dims - 1);

  points
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
  num_points = size(points,1);

  #creating the sparse vector mapping indices -> node ids
  ind2node = sparsevec( Int[sub2ind(max_dims, points[i,:]... ) for i=1:num_points ],
  	                    1:num_points,
  	                    prod(max_dims),
  	                    (x,y) -> x #ignore duplicates
  	                    )

  ind2node, max_dims
end


"""

    compute_DBF( point_cloud )

  Returns an array of DBF values for the point cloud. Currently creates
  a binary image, and runs bwd2 on it, though ideally we'd get rid of the
  need for an explicit bin_im

  TO OPTIMIZE
"""
function compute_DBF( point_cloud )

  bin_im = create_boundary_image( point_cloud );

  dbf_im = BWDists.bwd2( bin_im );

  dbf_vec = extract_dbf_values( dbf_im, point_cloud );

  dbf_vec;
end


"""

    create_boundary_image( point_cloud )

  Creates a boolean volume where the non-segment indices
  map to true, while the segment indices map to false.
"""
function create_boundary_image( point_cloud );

  max_dims = maximum( point_cloud, 1 );

  bin_im = ones(Bool, max_dims...);

  for p in 1:size( point_cloud, 1 )
    bin_im[ point_cloud[p,:]... ] = false;
  end

  bin_im;
end


"""

    extract_dbf_values( dbf_image, point_cloud )

  Takes an array where rows indicate subscripts, and extracts the values
  within a volume at those subscripts (in row order)
"""
@generated function extract_dbf_values{N}( dbf_image::Array{Float64,N}, point_cloud )
  quote

  num_points = size( point_cloud, 1 );
  dbf_values = zeros(num_points);

  for p in 1:size(point_cloud, 1)
    dbf_values[p] = (@nref $N dbf_image i->point_cloud[p,i]);
  end

  dbf_values
  end#quote
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
function make_neighbor_graph{T}( points::Array{T,2}, ind2node=nothing, max_dims=nothing )

  if ind2node == nothing ind2node, max_dims = create_node_lookup(points) end

  #init
  num_points = size(points,1);
  g = LightGraphs.Graph(num_points);

  #26-connectivity neighborhood
  # weights computed by euc_dist to center voxel
  nhood = [[i,j,k] for i=1:3,j=1:3,k=1:3];
  map!( x-> x - [2,2,2], nhood );
  nhood_weights = map( x-> sqrt(sum( x.^2 )), nhood );

  #only adding weights for non-duplicate nodes
  _,_,non_duplicates = findnz( ind2node );

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

    alexs_penalty( weights, DBF )

  Returns a version of the edge weights, modified by the DBF-based
  teasar penalty:

  w = w * 5000 .* (1 - DBF/maxDBF).^(16)

  The factor of 5000 shouldn't make a difference, but I've
  left it here
"""
function alexs_penalty( weights, DBF, G, mult=true )

  # *1.01 ensures later quotient is >0
  M = maximum(DBF,1)*1.01

  p_v = (1 - DBF./M ).^(16);

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

    local_max_multiplicative_penalty( weights, DBF, G )

  Returns a version of the edge weights, modified by the DBF-based
  penalty

  w = w * (1 - DBFdest/DBFstar)

  Where DBFstar is the maximum DBF value of an outwards neighbor for
  each node
"""
function local_max_multiplicative_penalty( weights, DBF, G )

  ns = LightGraphs.vertices(G)
  DBF_star = Vector{Float64}( length(ns) );

  #this is ok since ns is always an int range 1:N
  for n in ns
    neighborhood = LightGraphs.out_neighbors(G, n);

    if length(neighborhood) > 0
      DBF_star[n] = maximum(DBF[neighborhood]);
    else
      DBF_star[n] = 1
    end
  end

  is, js = findn(weights)
  ws = nonzeros(weights)
  new_ws = Vector{Float64}(length(ws));

  for i in eachindex(is)
    DBF_dest = DBF[js[i]];
    DBF_s = DBF_star[is[i]];

    new_ws[i] = ws[i] * (1 - DBF_dest/DBF_s);
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
  scale_param::Int=6, const_param::Int=6 );

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
function nodes_within_radius( sub::Array{Int,2}, ind2node, r, max_dims::Vector );

  beginning = convert(Vector{Int}, ceil(max(sub[:] .- r,1)));
  ending    = convert(Vector{Int}, floor(min(sub[:] .+ r, max_dims)));
  ind::Int = 0;

  nodes = Set{Int}();
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
function consolidate_paths( path_list )

  nodes = Set{Int}();
  edges = Set{Tuple{Int,Int}}();

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

end#module
