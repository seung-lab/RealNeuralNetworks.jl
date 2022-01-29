var documenterSearchIndex = {"docs":
[{"location":"man/references/","page":"-","title":"-","text":"Sümbül U, Song S, McCulloch K, Becker M, Lin B, Sanes JR, Masland RH, Seung HS. A genetic and computational approach to structurally classify neuronal types. Nature communications. 2014 Mar 24;5:3512. link\nSchierwagen A, Villmann T, Alpár A, Gärtner U. Cluster analysis of cortical pyramidal neurons using som. InIAPR Workshop on Artificial Neural Networks in Pattern Recognition 2010 Apr 11 (pp. 120-130). Springer, Berlin, Heidelberg.\nCuntz H, Forstner F, Borst A, H\\äusser M. The TREES Toolbox—Probing the Basis of Axonal and Dendritic Segmenting. Neuroinformatics. 2011;1–6. \nSchierwagen A. Neuronal morphology: Shape characteristics and models. Neurophysiology. 2008;40(4):310–315. \nUylings HB., van Pelt J. Measures for quantifying dendritic arborizations. Network: Computation in Neural Systems. 2002;13(3):397–414. a good review paper\nWanner AA, Genoud C, Masudi T, Siksou L, Friedrich RW. Dense EM-based reconstruction of the interglomerular projectome in the zebrafish olfactory bulb. Nature neuroscience. 2016 Jun 1;19(6):816-25. also have clustering methods\nSato, M., I. Bitter, M. A. Bender, A. E. Kaufman, and M. Nakajima. “TEASAR: Tree-Structure Extraction Algorithm for Accurate and Robust Skeletons.” In Proceedings the Eighth Pacific Conference on Computer Graphics and Applications, 281–449, 2000. doi:10.1109/PCCGA.2000.883951.","category":"page"},{"location":"lib/internals/PointArrays/#PointArrays","page":"PointArrays","title":"PointArrays","text":"","category":"section"},{"location":"lib/internals/PointArrays/","page":"PointArrays","title":"PointArrays","text":"Modules = [RealNeuralNetworks.NodeNets.PointArrays]","category":"page"},{"location":"lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.add_offset!-Union{Tuple{T}, Tuple{Matrix{T}, Tuple{T, T, T}}} where T","page":"PointArrays","title":"RealNeuralNetworks.NodeNets.PointArrays.add_offset!","text":"add offset to points\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.from_binary_image-Tuple{Union{BitArray{3}, Array{Bool, 3}}}","page":"PointArrays","title":"RealNeuralNetworks.NodeNets.PointArrays.from_binary_image","text":"parameter:     bin_im: binary array. object voxels are false, non-object voxels are true!\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.from_seg-Union{Tuple{Array{T, 3}}, Tuple{T}} where T","page":"PointArrays","title":"RealNeuralNetworks.NodeNets.PointArrays.from_seg","text":"find points inside an object from a segmentation array. \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.get_boundary_point_indexes-Union{Tuple{TSeg}, Tuple{T}, Tuple{Matrix{T}, Array{TSeg, 3}}, Tuple{Matrix{T}, Array{TSeg, 3}, TSeg}} where {T, TSeg}","page":"PointArrays","title":"RealNeuralNetworks.NodeNets.PointArrays.get_boundary_point_indexes","text":"find out the boundary voxels and represent them as indexes in the point array\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.merge-Union{Tuple{T}, Tuple{Matrix{T}, Matrix{T}}} where T","page":"PointArrays","title":"RealNeuralNetworks.NodeNets.PointArrays.merge","text":"merge(self::Array{T,3}, other::Array{T,3})\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/#Internal-Documentation","page":"Internals","title":"Internal Documentation","text":"","category":"section"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"This page lists all the documented internals of the RealNeuralNetworks module and submodules.","category":"page"},{"location":"lib/internals/#Contents","page":"Internals","title":"Contents","text":"","category":"section"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]","category":"page"},{"location":"lib/internals/#Index","page":"Internals","title":"Index","text":"","category":"section"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"A list of all internal documentation sorted by module. ","category":"page"},{"location":"lib/internals/","page":"Internals","title":"Internals","text":"Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]","category":"page"},{"location":"man/contributing/","page":"-","title":"-","text":"Please fork and create pull requests in github.","category":"page"},{"location":"man/getting_started/","page":"Getting Started","title":"Getting Started","text":"using RealNeuralNetworks.NodeNets\nusing RealNeuralNetworks.Neurons\nusing RealNeuralNetworks.SWCs\n\n# skeletonization\nnodeNet = NodeNet(seg::Array{UInt32,3}; obj_id = convert(UInt32,77605))\nneuron = Neuron( nodeNet )\nswc = SWC(neuron)\nSWCs.save(swc, tempname()*\".swc\")","category":"page"},{"location":"lib/internals/NodeNets/#NodeNets","page":"NodeNets","title":"NodeNets","text":"","category":"section"},{"location":"lib/internals/NodeNets/","page":"NodeNets","title":"NodeNets","text":"Modules = [\n    RealNeuralNetworks.NodeNets,\n    RealNeuralNetworks.NodeNets.DBFs\n]","category":"page"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.NodeNet-Tuple{Union{Array{Bool, 3}, BitArray}}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.NodeNet","text":"NodeNet(bin_im)\n\nParameters:     bin_im: binary mask. the object voxel should be false, non-object voxel should be true Return:     nodeNet object\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.NodeNet-Union{Tuple{Array{T, 3}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.NodeNet","text":"NodeNet( seg, obj_id; penalty_fn=alexs_penalty)\n\nPerform the teasar algorithm on the passed binary array.\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.NodeNet-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.NodeNet","text":"NodeNet( points; penalty_fn = alexs_penalty )\n\nPerform the teasar algorithm on the passed Nxd array of points\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.alexs_penalty","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.alexs_penalty","text":"alexs_penalty( weights, dbf )\n\nReturns a version of the edge weights, modified by the DBF-based   teasar penalty:\n\nw = w * 5000 .* (1 - dbf/maxDBF).^(16)\n\nThe factor of 5000 shouldn't make a difference, but I've   left it here\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.consolidate_paths-Tuple{Vector}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.consolidate_paths","text":"consolidate_paths( path_list )\n\nExtracts the unique nodes and edges among all paths in the list\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.create_node_lookup-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.create_node_lookup","text":"create_node_lookup( points )\n\nAbstractly represent the points as a volume by means of   linear indexes into a sparse vector.\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.distill!-Union{Tuple{T}, Tuple{Matrix{T}, Vector, Vector}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.distill!","text":"Parameters:\n\npoint_array: a Nx3 array recording all the voxel coordinates inside the object.\npath_nodes: the indexes of nodeNet voxels in the point_array \npath_edges: the index pairs of nodeNet in the point_array\n\nIn the pathnodes, the nodes were encoded as the index of the pointarray. Since we do  not need pointarray anymore, which could be pretty big, we'll only reserve the nodeNet  coordinates and let the gc release the pointarray.  the same applys to path_edges \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.find_closest_node_id-Union{Tuple{T}, Tuple{N}, Tuple{RealNeuralNetworks.NodeNets.NodeNet{T}, Tuple{Vararg{T, N}}}} where {N, T}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.find_closest_node_id","text":"find_closest_node_id(self::NodeNet{T}, point::NTuple{3,T}) where T\n\nlook for the id of the closest node\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.find_new_root_node_V1-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.find_new_root_node_V1","text":"find_new_root_node( points )\n\nExtracts the point with the lowest linear index from the   Array of passed points \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.find_new_root_node_id-Tuple{Vector{Float32}, DataStructures.IntSet}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.find_new_root_node_id","text":"find_new_root_node( dbf )\n\nExtracts the point with the largest dbf from the    Array of passed points \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_connectivity_matrix-Tuple{Vector}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.get_connectivity_matrix","text":"get_connectivity_matrix(edges::Vector)\n\nconstruct sparse connectivity matrix accordint to the edges \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_edge_num-Tuple{RealNeuralNetworks.NodeNets.NodeNet}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.get_edge_num","text":"the connectivity matrix is symmetric, so the connection is undirected\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_neuroglancer_precomputed-Tuple{RealNeuralNetworks.NodeNets.NodeNet}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.get_neuroglancer_precomputed","text":"get binary buffer formatted as neuroglancer nodeNet.\n\nBinary format\n\nUInt32: number of vertex\nUInt32: number of edges\nArray{Float32,2}: Nx3 array, xyz coordinates of vertex\nArray{UInt32,2}: Mx2 arrray, node index pair of edges\n\nreference:  https://github.com/seung-lab/neuroglancer/wiki/Skeletons\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_node_array-Tuple{RealNeuralNetworks.NodeNets.NodeNet}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.get_node_array","text":"get_node_array(self::NodeNet)\n\nReturn:     nodeArray::Array{T,2}: size is (nd, np), nd is the dimention of each node,                  np is the number of nodes\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_num_segment_point-Tuple{RealNeuralNetworks.NodeNets.NodeNet}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.get_num_segment_point","text":"get_segment_point_num(self::NodeNet)\n\nget number of branching points \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_num_segmentes-Tuple{RealNeuralNetworks.NodeNets.NodeNet}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.get_num_segmentes","text":"assume that the graph is acyclic, no loop.\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_sholl_number-Tuple{RealNeuralNetworks.NodeNets.NodeNet, AbstractFloat}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.get_sholl_number","text":"get_sholl_number(self::NodeNet, radius::AbstractFloat)\n\nget the number of points which is in neurite and incounters with a sphere centered on root node \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.local_max_multiplicative_penalty-Tuple{Any, Any, Any}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.local_max_multiplicative_penalty","text":"local_max_multiplicative_penalty( weights, dbf, G )\n\nReturns a version of the edge weights, modified by the DBF-based   penalty\n\nw = w * (1 - DBFdest/DBFstar)\n\nWhere DBFstar is the maximum DBF value of an outwards neighbor for   each node\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.make_neighbor_graph-Union{Tuple{Matrix{T}}, Tuple{T}, Tuple{Matrix{T}, Any}, Tuple{Matrix{T}, Any, Any}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.make_neighbor_graph","text":"make_neighbor_graph( points )\n\nCreates a LightGraphs graph from a point cloud, such that   each point has a node (indexed in the same order as the rows),   and each node is connected to its neighbors (by 26-connectivity).\n\nAlso returns a sparse matrix of weights for these edges equal to the   euclidean distance between the indices of each point. These can   be weighted or modified easily upon return.\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.nodes_within_radius-Union{Tuple{T}, Tuple{Vector{T}, Any, Any, Vector}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.nodes_within_radius","text":"nodes_within_radius( sub, sub2node, r, max_dims )\n\nIdentifies the node indices within the subscript\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.remove_path_from_rns!","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.remove_path_from_rns!","text":"remove_path_from_rns!( reachableNodeList, path, points, sub2node, dbf, max_dims, scale_param, const_param)\n\nIdentifies the nodes to remove from the reachable nodes within the graph.   Probably the ugliest function here, and should be optimized later\n\nTO OPTIMIZE\n\n\n\n\n\n","category":"function"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.save-Tuple{RealNeuralNetworks.NodeNets.NodeNet, Integer, AbstractDict}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.save","text":"save(self::NodeNet, cellId::UInt32, d_json::Associative, d_bin::Associative)\n\nsave nodeNet in google cloud storage for neuroglancer visualization the format is the same with meshes\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.save_edges-Tuple{RealNeuralNetworks.NodeNets.NodeNet, String}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.save_edges","text":"save binary file of point pairs used in neuroglancer python interface to visualize the nodeNet \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.translate_to_origin!-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.translate_to_origin!","text":"translate_to_origin!( points::Matrix{T} )\n\nNormalize the point dimensions by subtracting the min   across each dimension. This step isn't extremely necessary,   but might be useful for compatibility with the MATLAB code.   record the offset and add it back after building the nodeNet\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.compute_DBF","text":"compute_DBF( pointCloud )\n\nReturns an array of DBF values for the point cloud. Currently creates   a binary image, and runs bwd2 on it, though ideally we'd get rid of the   need for an explicit bin_im\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{T}, Tuple{Array{T, 3}, T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.compute_DBF","text":"use segmentation to get binary image to save memory usage\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{T}, Tuple{Matrix{T}, Union{Array{Bool, 3}, BitArray}}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.compute_DBF","text":"compute_DBF( bin_im )\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{T}, Tuple{Matrix{T}, Vector}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.compute_DBF","text":"compute Distance from Boundary Field (DBF) based on point cloud and the boundary points\n\nWARN: this function do not work correctly!\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.create_binary_image-Union{Tuple{Array{T, 3}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.create_binary_image","text":"create_binary_image( seg, obj_id )\n\nCreates a boolean volume where the non-segment indices map to true, while the segment indices map to false \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.create_binary_image-Union{Tuple{Matrix{T}}, Tuple{T}} where T","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.create_binary_image","text":"create_binary_image( pointCloud )\n\nCreates a boolean volume where the non-segment indices   map to true, while the segment indices map to false.\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.distance_transform-Union{Tuple{AbstractArray{T, N}}, Tuple{N}, Tuple{T}, Tuple{AbstractArray{T, N}, Vector{Float32}}} where {T, N}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.distance_transform","text":"distance_transform( d::AbstractArray{T,N}, voxelSize::Vector{Float32}=ones(Float32, N) )\n\nReturns a euclidean distance transformation of the mask provided by d. The return   value will be a volume of the same size as d where the value at each index corresponds   to the distance between that location and the nearest location for which d > 0.\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.extract_dbf_values-Union{Tuple{N}, Tuple{Array{Float32, N}, Any}} where N","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.extract_dbf_values","text":"extract_dbf_values( dbf_image, pointCloud )\n\nTakes an array where rows indicate subscripts, and extracts the values   within a volume at those subscripts (in row order)\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.fill_f0!-Union{Tuple{N}, Tuple{T}, Tuple{Array{Float32, N}, AbstractArray{T, N}}} where {T, N}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.fill_f0!","text":"Fills an n-dimensional volume with initial states for edt transformation, inf for non-feature voxels, and 0 for feature voxels\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.fv_isfurther-Tuple{Float32, Int64, Float32, Int64, Int64}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.fv_isfurther","text":"Getting too tired to document these next few, but will be worth it if it works\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.remove_euclidean_distance_transform-Tuple{Float32, Float32, Float32, Int64, Int64, Int64}","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.remove_euclidean_distance_transform","text":"remove_euclidean_distance_transform\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.row_voronoi_edt!-Union{Tuple{N}, Tuple{Array{Float32, N}, Tuple, Vector{Float32}, Vector{Int64}, Float32}} where N","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.row_voronoi_edt!","text":"Performs the edt over a specific row in the volume, following the first dimension\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.vol_voronoi_edt!-Union{Tuple{N}, Tuple{Array{Float32, N}, Float32}} where N","page":"NodeNets","title":"RealNeuralNetworks.NodeNets.DBFs.vol_voronoi_edt!","text":"Performs the edt transformation along the first dimension of the N-dimensional volume\n\n\n\n\n\n","category":"method"},{"location":"lib/public/#Public-Documentation","page":"Public","title":"Public Documentation","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"Documentation for RealNeuralNetworks.jl's public interface. ","category":"page"},{"location":"lib/public/","page":"Public","title":"Public","text":"See Internal Documentation for internal package docs covering all submodules. ","category":"page"},{"location":"lib/public/#Contents","page":"Public","title":"Contents","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#Index","page":"Public","title":"Index","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"Pages = [\"public.md\"]","category":"page"},{"location":"lib/public/#Public-Interface","page":"Public","title":"Public Interface","text":"","category":"section"},{"location":"lib/public/","page":"Public","title":"Public","text":"RealNeuralNetworks\nNeurons\nSegments\nSynapses\nSWCs\nNeuralNets\nSynapseTables","category":"page"},{"location":"#RealNeuralNetworks.jl","page":"Home","title":"RealNeuralNetworks.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"A Skeletonization, morphological and connectivity analysis tool. ","category":"page"},{"location":"#Features","page":"Home","title":"Features","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Skeletonization using TEASAR algorithm. \nSkeleton IO with format of SWC and neuroglancer precomputed. \nA lot of morphological analysis functions. \nSynaptic distribution analysis in/between neurons.\nNeuron connectivity analysis.","category":"page"},{"location":"","page":"Home","title":"Home","text":"The Getting Started provides a tutorial explaning how to get started.","category":"page"},{"location":"","page":"Home","title":"Home","text":"Some examples of packages could be found on the Examples page. ","category":"page"},{"location":"","page":"Home","title":"Home","text":"See the Index for the complete list of functions and types. ","category":"page"},{"location":"#Manual-Outline","page":"Home","title":"Manual Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"man/getting_started.md\",\n    \"man/examples.md\",\n    <!-- \"man/contributing.md\", -->\n    \"man/references.md\"\n]\nDepth=1","category":"page"},{"location":"#Library-Outline","page":"Home","title":"Library Outline","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"lib/public.md\",\n    \"lib/internals.md\"\n]\n\n###[Index](@id main-index)\n","category":"page"},{"location":"","page":"Home","title":"Home","text":"@index  Pages = [\"lib/public.md\"] ```","category":"page"},{"location":"man/examples/#Examples","page":"Examples","title":"Examples","text":"","category":"section"},{"location":"man/examples/","page":"Examples","title":"Examples","text":"All the examples are hosted in another GitHub Repo using Jupyter Notebooks","category":"page"},{"location":"lib/internals/Manifests/#Manifests","page":"Manifests","title":"Manifests","text":"","category":"section"},{"location":"lib/internals/Manifests/","page":"Manifests","title":"Manifests","text":"Modules = [\n    RealNeuralNetworks.Manifests,\n]","category":"page"},{"location":"lib/internals/Manifests/#RealNeuralNetworks.Manifests.Manifest-Tuple{AbstractString, AbstractString, AbstractString, Integer}","page":"Manifests","title":"RealNeuralNetworks.Manifests.Manifest","text":"Parameters:     dir \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/Manifests/#RealNeuralNetworks.Manifests.Manifest-Union{Tuple{T}, Tuple{D}, Tuple{Vector, BigArrays.BigArray{D, T}}} where {D, T}","page":"Manifests","title":"RealNeuralNetworks.Manifests.Manifest","text":"example: [\"770048087:0:2968-34801776-228816912-17424\"]\n\n\n\n\n\n","category":"method"},{"location":"lib/internals/Manifests/#RealNeuralNetworks.Manifests.get_voxel_offset-Tuple{RealNeuralNetworks.Manifests.Manifest}","page":"Manifests","title":"RealNeuralNetworks.Manifests.get_voxel_offset","text":"the voxel offset in the neuroglancer precomputed info file \n\n\n\n\n\n","category":"method"},{"location":"lib/internals/Manifests/#RealNeuralNetworks.Manifests.trace-Tuple{RealNeuralNetworks.Manifests.Manifest}","page":"Manifests","title":"RealNeuralNetworks.Manifests.trace","text":"iterate the chunks containing the neuron with specified neuronId build point cloud and dbf when iterating the chunks \n\n\n\n\n\n","category":"method"}]
}
