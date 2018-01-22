module NodeNets
include("DBFs.jl"); using .DBFs;
include("PointArrays.jl"); using .PointArrays;
using ..RealNeuralNetworks.SWCs
using Base.Cartesian

export NodeNet 

using LightGraphs
using MetaGraphs 
import BigArrays

const ZERO_UINT32 = convert(UInt32, 0)
const ONE_UINT32  = convert(UInt32, 1)
const OFFSET = (ZERO_UINT32, ZERO_UINT32, ZERO_UINT32)

# rescale the skeleton
const EXPANSION = (ONE_UINT32, ONE_UINT32, ONE_UINT32)

# control the removing points around path based on DBF
const REMOVE_PATH_SCALE = 3 
const REMOVE_PATH_CONST = 4

const NodeNet = MetaGraph

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
function NodeNet(point_array::Array{T,2}, path_nodes::Vector{Int}, 
                         path_edges::Vector, node_radii::Vector) where T 
    numNodes = length(path_nodes) 
    numEdges = length(path_edges)
    nodeNet = MetaGraph( numNodes )
    # map the old path node id to new id 
    id_map = Dict{T, T}()
    sizehint!(id_map, numNodes)
    # build new nodes and the id map
    for i in 1:numNodes 
        coordinate = (point_array[path_nodes[i], :]...)
        set_props!(nodeNet, i, Dict(:coordinate => Float32.(coordinate), 
                                    :radius     => Float32(node_radii[i])))
        id_map[path_nodes[i]] = i
    end 

    # rebuild the edges
    for i in 1:numEdges
        edge = (id_map[path_edges[i][1]], id_map[path_edges[i][2]])
        @assert edge[1] != edge[2]
        add_edge!(nodeNet, edge)
    end
    nodeNet
end 

"""
    NodeNet( points; penalty_fn = alexs_penalty )

  Perform the teasar algorithm on the passed Nxd array of points
"""
function NodeNet( points::Array{T,2}; dbf=DBFs.compute_DBF(points),
                         penalty_fn::Function = alexs_penalty,
                         expansion::NTuple{3, UInt32} = EXPANSION) where T

    println("total number of points: $(size(points,1))")
    points, bbox_offset = translate_to_origin!( points );
    ind2node, max_dims = create_node_lookup( points );
    sub2node = x -> ind2node[ sub2ind(max_dims, x[1],x[2],x[3]) ];   

    println("making graph (2 parts)");
    @time euclideanWeightedGraph = make_neighbor_graph( points, ind2node, max_dims;)
    println("build dbf weights from penalty function ...")
    @time dbfWeightedGraph = penalty_fn( euclideanWeightedGraph, dbf )
 
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
        println("root: ", root_node);
 
        # do a graph traversal to find farthest node and
        # find reachable nodes
        dsp_euclidean = LightGraphs.dijkstra_shortest_paths(euclideanWeightedGraph, 
                                                            root_node);
        # and another to find the min DBF-weighted paths to all nodes
        dsp_dbf = LightGraphs.dijkstra_shortest_paths(dbfWeightedGraph, root_node);
    
        # remove reachable nodes from the disconnected ones
        reachable_nodes = find(.!(isinf.(dsp_euclidean.dists)));
        # another necessary precaution for duplicate nodes
        reachable_nodes = intersect(reachable_nodes, disconnected_nodes);
        setdiff!(disconnected_nodes, reachable_nodes);
        empty!(inspected_nodes);
    
        while length(reachable_nodes) > 0
            # find the node farthest away from the root by euc distance
            _, farthest_index = findmax( dsp_euclidean.dists[[reachable_nodes...]] );
            farthest_node = reachable_nodes[farthest_index];
            push!(destinations, farthest_node);
            println("dest: ", farthest_node)
    
            if farthest_node == root_node break end # this can happen apparently
    
            new_path = LightGraphs.enumerate_paths( dsp_dbf, farthest_node );
    
            push!(paths, new_path)
            #can't do this in-place with arrays
            #this fn call is getting ridiculous
            @time reachable_nodes = remove_path_from_rns( reachable_nodes, new_path,
                                                    points, ind2node,
                                                    dbf, [max_dims...],
                                                    inspected_nodes );
        end #while reachable nodes from root
    end #while disconnected nodes
  
    println("Consolidating Paths...")
    path_nodes, path_edges = consolidate_paths( paths );
    node_radii = dbf[path_nodes];
    
    println("build new NodeNet...")
    nodeNet = NodeNet(points, path_nodes, path_edges, node_radii)
    add_offset!(nodeNet, bbox_offset)
    return nodeNet
end

function NodeNet(swc::SWC)
    numNodes = length(swc)
    nodeNet = NodeNet(numNodes) 
    for (index, point) in enumerate(swc)
        set_props!(nodeNet, index, 
                        Dict(:coordinate => (point.x, point.y, point.z), 
                             :radius     => point.radius))
        if point.parent != -1
            add_edge!(nodeNet, (index, point.parent))
        end 
    end 
    nodeNet
end 

function SWCs.SWC(nodeNet::NodeNet)
    swc = SWC()
    sizehint!(swc, NodeNets.get_node_num(nodeNet))
    
    for v in vertices(nodeNet)
        coordinate = get_prop(nodeNet, v, :coordinate)
        radius = get_prop(nodeNet, v, :radius)
        point = SWCs.PointObj(0, coordinate..., radius, -1)
        push!(swc, point)
    end
    # assign parents according to edge 
    for edge in edges(nodeNet)
        swc[dst(edge)].parent = src(edge)
    end  
    swc
end

##################### properties ###############################
function get_node_list(self::NodeNet) 
    N = get_node_num(self)
    nodeList = Vector{NTuple{3, Float32}}()
    sizehint!(nodeList, N)
    for v in vertices(self)
        coordinate = get_prop(self, v, :coordinate)
        push!(nodeList, coordinate)
    end 
    nodeList 
end 
@inline function get_connectivity_matrix(self::NodeNet) adjacency_matrix(self) end
@inline function get_radii(self::NodeNet) 
    map(v->get_prop(self,v,:radius),  vertices(self)) 
end 
@inline function get_node_num(self::NodeNet) length(vertices(self)) end
@inline function get_edge_num(self::NodeNet) length(edges(self)) end
@inline function get_edges(self::NodeNet) edges(self) end  

function get_num_of_branching_points(self::NodeNet)
    conn = adjacency_matrix(self)
    num_branching_point = 0
    for col in size(conn, 2)
        if nnz(conn[:,col]) > 2
            num_branching_point += 1
        end
    end 
    num_branching_point 
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
        write(buffer, UInt32( src(edge) - 1 ))
        write(buffer, UInt32( dst(edge) - 1 ))
    end
    bin = Vector{UInt8}(take!(buffer))
    close(buffer)
    return bin 
end 

##################### manipulate ############################
@inline function add_offset!(self::NodeNet, offset::Tuple )
    @assert length(offset) == 3
    for v in vertices(self)
        coordinate = get_prop(self, v, :coordinate)
        coordinate = map(+, coordinate, Float32.(offset))
        set_prop!(self, v, :coordinate, coordinate)
    end
    nothing
end

@inline function stretch_coordinates!(self::NodeNet, mip::Integer)
    expansion = (2.0^mip, 2.0^mip, 1.0) 
    stretch_coordinates!(self, expansion)
end 
@inline function stretch_coordinates!(self::NodeNet, expansion::Tuple)
    @assert length(expansion) == 3
    for v in vertices(self)
        coordinate = get_prop(self, v, :coordinate)
        coordinate = map(*, coordinate, Float32.(expansion))
        set_prop!(self, v, :coordinate, coordinate)
    end       
    nothing
end 

#---------------------------------------------------------------


"""

    translate_to_origin( points )

  Normalize the point dimensions by subtracting the min
  across each dimension. This step isn't extremely necessary,
  but might be useful for compatibility with the MATLAB code.
  record the offset and add it back after building the nodeNet
"""
@inline function translate_to_origin!( points::Array{T,2} ) where T
    offset = minimum( points, 1 ) - T(1) ;
    # transform to 1d vector
    offset = (vec( offset )...)
    @assert length(offset) == 3
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
function create_node_lookup( points::Array{T,2} ) where T
  # need the max in order to define the bounds of that volume
  # tuple allows the point to be passed to fns as dimensions
  max_dims = ( round.(Int, maximum( points, 1 ))... );
  num_points = size(points,1);
  #creating the sparse vector mapping indices -> node ids
  ind2node = sparsevec( [sub2ind(max_dims, Vector{Int}(points[i,:])... ) 
                                                        for i=1:num_points ],
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
function make_neighbor_graph( points::Array{T,2}, ind2node=nothing, max_dims=nothing;
                             expansion::NTuple{3, UInt32} = EXPANSION ) where T

    if ind2node == nothing ind2node, max_dims = create_node_lookup(points) end
  
    #init
    numPoints = size(points,1);
    g = MetaGraph(numPoints, 0.0);
  
    #26-connectivity neighborhood
    # weights computed by euc_dist to center voxel
    nhood = [[i,j,k] for i=1:3,j=1:3,k=1:3];
    map!( x-> (x .- [2,2,2]).*[expansion ...] , nhood, nhood)
    const nhood_weights = map( norm, nhood );
  
    #only adding weights for non-duplicate nodes
    _,non_duplicates = findnz( ind2node );
  
    #init - descriptors of the edges for graph, weights later
    is = Vector{Int}();
    js = Vector{Int}();
  
    dest_sub = zeros(Int,3)
  
    for i in eachindex(nhood)
    	# only need to explicitly tell LightGraphs to add the edges in
    	# one direction (undirected), though doing this properly can be tricky
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
            dest_ind = sub2ind( max_dims, dest_sub... );
            dest_node = ind2node[dest_ind];
  
            if dest_node != 0 #valid edge
                edge = LightGraphs.Edge(p, dest_node)
                LightGraphs.add_edge!(g, edge)
                set_prop!(g, edge, :weight, nhood_weights[i])
                push!(is,p); push!(js,dest_node); 
            end
  
        end #for p
    end #for i in nhood
  
    g 
end


"""

    alexs_penalty( weights, dbf )

  Returns a version of the edge weights, modified by the DBF-based
  teasar penalty:

  w = w * 5000 .* (1 - dbf/maxDBF).^(16)

  The factor of 5000 shouldn't make a difference, but I've
  left it here
"""
function alexs_penalty( G::MetaGraph, dbf::Vector, mult::Bool=true )

    # *1.01 ensures later quotient is >0
    const M = maximum(dbf,1)*1.01

    const p_v = (1 - dbf./M ).^(16);

    #DBF penalty is defined by the destination node, only need js
    weights = LightGraphs.weights(G)
    dbfWeightedGraph = copy(G)
    for edge in edges(G)
        dst = LightGraphs.dst(edge)
        weight = get_prop(G, edge, :weight)
        set_prop!(dbfWeightedGraph, edge, :weight, weight*p_v[dst])
    end 
    dbfWeightedGraph
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

    remove_path_from_rns( reachable_nodes, path, points, sub2node, dbf, max_dims, scale_param, const_param)

  Identifies the nodes to remove from the reachable nodes within the graph.
  Probably the ugliest function here, and should be optimized later

  TO OPTIMIZE
"""
function remove_path_from_rns( reachable_nodes::Vector, path::Vector,
  points, ind2node, dbf, max_dims, inspected_nodes,
  scale_param::Int=REMOVE_PATH_SCALE, const_param::Int=REMOVE_PATH_CONST );

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
function nodes_within_radius( sub::Array{T,1}, ind2node, r, max_dims::Vector ) where T;

  beginning = convert(Vector{Int}, ceil.(max.(sub[:] .- r,1)));
  ending    = convert(Vector{Int}, floor.(min.(sub[:] .+ r, max_dims)));
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

end # module
