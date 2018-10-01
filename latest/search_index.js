var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#RealNeuralNetworks.jl-1",
    "page": "Home",
    "title": "RealNeuralNetworks.jl",
    "category": "section",
    "text": "A Skeletonization, morphological and connectivity analysis tool. "
},

{
    "location": "#Features-1",
    "page": "Home",
    "title": "Features",
    "category": "section",
    "text": "Skeletonization using TEASAR algorithm. \nSkeleton IO with format of SWC and neuroglancer precomputed. \nA lot of morphological analysis functions. \nSynaptic distribution analysis in/between neurons.\nNeuron connectivity analysis.The Getting Started provides a tutorial explaning how to get started.Some examples of packages could be found on the Examples page. See the Index for the complete list of functions and types. "
},

{
    "location": "#Manual-Outline-1",
    "page": "Home",
    "title": "Manual Outline",
    "category": "section",
    "text": "Pages = [\n    \"man/getting_started.md\",\n    \"man/examples.md\",\n    \"man/contributing.md\",\n    \"man/references.md\"\n]\nDepth=1"
},

{
    "location": "#Library-Outline-1",
    "page": "Home",
    "title": "Library Outline",
    "category": "section",
    "text": "Pages = [\n    \"lib/public.md\",\n    \"lib/internals.md\"\n]\n\n###[Index](@id main-index)\n@index  Pages = [\"lib/public.md\"] ```"
},

{
    "location": "man/getting_started/#",
    "page": "Getting Started",
    "title": "Getting Started",
    "category": "page",
    "text": "using RealNeuralNetworks.NodeNets\nusing RealNeuralNetworks.Neurons\nusing RealNeuralNetworks.SWCs\n\n# skeletonization\nnodeNet = NodeNet(seg::Array{UInt32,3}; obj_id = convert(UInt32,77605))\nneuron = Neuron( nodeNet )\nswc = SWC(neuron)\nSWCs.save(swc, tempname()*\".swc\")"
},

{
    "location": "man/examples/#",
    "page": "Examples",
    "title": "Examples",
    "category": "page",
    "text": ""
},

{
    "location": "man/examples/#Examples-1",
    "page": "Examples",
    "title": "Examples",
    "category": "section",
    "text": "show some usage examples. "
},

{
    "location": "man/references/#",
    "page": "-",
    "title": "-",
    "category": "page",
    "text": "Sümbül U, Song S, McCulloch K, Becker M, Lin B, Sanes JR, Masland RH, Seung HS. A genetic and computational approach to structurally classify neuronal types. Nature communications. 2014 Mar 24;5:3512. link\nSchierwagen A, Villmann T, Alpár A, Gärtner U. Cluster analysis of cortical pyramidal neurons using som. InIAPR Workshop on Artificial Neural Networks in Pattern Recognition 2010 Apr 11 (pp. 120-130). Springer, Berlin, Heidelberg.\nCuntz H, Forstner F, Borst A, H\\äusser M. The TREES Toolbox—Probing the Basis of Axonal and Dendritic Segmenting. Neuroinformatics. 2011;1–6. \nSchierwagen A. Neuronal morphology: Shape characteristics and models. Neurophysiology. 2008;40(4):310–315. \nUylings HB., van Pelt J. Measures for quantifying dendritic arborizations. Network: Computation in Neural Systems. 2002;13(3):397–414. a good review paper\nWanner AA, Genoud C, Masudi T, Siksou L, Friedrich RW. Dense EM-based reconstruction of the interglomerular projectome in the zebrafish olfactory bulb. Nature neuroscience. 2016 Jun 1;19(6):816-25. also have clustering methods\nSato, M., I. Bitter, M. A. Bender, A. E. Kaufman, and M. Nakajima. “TEASAR: Tree-Structure Extraction Algorithm for Accurate and Robust Skeletons.” In Proceedings the Eighth Pacific Conference on Computer Graphics and Applications, 281–449, 2000. doi:10.1109/PCCGA.2000.883951."
},

{
    "location": "lib/public/#",
    "page": "Public",
    "title": "Public",
    "category": "page",
    "text": ""
},

{
    "location": "lib/public/#Public-Documentation-1",
    "page": "Public",
    "title": "Public Documentation",
    "category": "section",
    "text": "Documentation for RealNeuralNetworks.jl\'s public interface. See Internal Documentation for internal package docs covering all submodules. "
},

{
    "location": "lib/public/#Contents-1",
    "page": "Public",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public/#Index-1",
    "page": "Public",
    "title": "Index",
    "category": "section",
    "text": "Pages = [\"public.md\"]"
},

{
    "location": "lib/public/#Public-Interface-1",
    "page": "Public",
    "title": "Public Interface",
    "category": "section",
    "text": "RealNeuralNetworks\nNeurons\nSegments\nSynapses\nSWCs\nNeuralNets\nSynapseTables"
},

{
    "location": "lib/internals/#",
    "page": "Internals",
    "title": "Internals",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals/#Internal-Documentation-1",
    "page": "Internals",
    "title": "Internal Documentation",
    "category": "section",
    "text": "This page lists all the documented internals of the RealNeuralNetworks module and submodules."
},

{
    "location": "lib/internals/#Contents-1",
    "page": "Internals",
    "title": "Contents",
    "category": "section",
    "text": "Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]"
},

{
    "location": "lib/internals/#Index-1",
    "page": "Internals",
    "title": "Index",
    "category": "section",
    "text": "A list of all internal documentation sorted by module. Pages = [joinpath(\"internals\", f) for f in readdir(\"internals\")]"
},

{
    "location": "lib/internals/Manifests/#",
    "page": "Manifests",
    "title": "Manifests",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals/Manifests/#RealNeuralNetworks.Manifests.Manifest-Tuple{AbstractString,AbstractString,AbstractString}",
    "page": "Manifests",
    "title": "RealNeuralNetworks.Manifests.Manifest",
    "category": "method",
    "text": "Parameters:     dir \n\n\n\n\n\n"
},

{
    "location": "lib/internals/Manifests/#RealNeuralNetworks.Manifests.Manifest-Tuple{Dict{Symbol,Any},BigArrays.AbstractBigArray}",
    "page": "Manifests",
    "title": "RealNeuralNetworks.Manifests.Manifest",
    "category": "method",
    "text": "example: {\"fragments\": [\"770048087:0:2968-34801776-228816912-17424\"]}\n\n\n\n\n\n"
},

{
    "location": "lib/internals/Manifests/#RealNeuralNetworks.Manifests.Manifest-Union{Tuple{C}, Tuple{N}, Tuple{T}, Tuple{D}, Tuple{Array{T,1} where T,BigArray{D,T,N,C}}} where C where N where T where D",
    "page": "Manifests",
    "title": "RealNeuralNetworks.Manifests.Manifest",
    "category": "method",
    "text": "example: [\"770048087:0:2968-34801776-228816912-17424\"]\n\n\n\n\n\n"
},

{
    "location": "lib/internals/Manifests/#RealNeuralNetworks.Manifests.get_voxel_offset-Tuple{RealNeuralNetworks.Manifests.Manifest}",
    "page": "Manifests",
    "title": "RealNeuralNetworks.Manifests.get_voxel_offset",
    "category": "method",
    "text": "the voxel offset in the neuroglancer precomputed info file \n\n\n\n\n\n"
},

{
    "location": "lib/internals/Manifests/#RealNeuralNetworks.Manifests.trace-Tuple{RealNeuralNetworks.Manifests.Manifest,Any}",
    "page": "Manifests",
    "title": "RealNeuralNetworks.Manifests.trace",
    "category": "method",
    "text": "iterate the chunks containing the neuron with specified cellId build point cloud and dbf when iterating the chunks \n\n\n\n\n\n"
},

{
    "location": "lib/internals/Manifests/#Manifests-1",
    "page": "Manifests",
    "title": "Manifests",
    "category": "section",
    "text": "Modules = [\n    RealNeuralNetworks.Manifests,\n]"
},

{
    "location": "lib/internals/NodeNets/#",
    "page": "NodeNets",
    "title": "NodeNets",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.NodeNet-Tuple{Union{Array{Bool,3}, BitArray}}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.NodeNet",
    "category": "method",
    "text": "NodeNet(bin_im)\n\nParameters:     bin_im: binary mask. the object voxel should be false, non-object voxel should be true Return:     nodeNet object\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.NodeNet-Union{Tuple{Array{T,2}}, Tuple{T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.NodeNet",
    "category": "method",
    "text": "NodeNet( points; penalty_fn = alexs_penalty )\n\nPerform the teasar algorithm on the passed Nxd array of points\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.NodeNet-Union{Tuple{Array{T,3}}, Tuple{T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.NodeNet",
    "category": "method",
    "text": "NodeNet( seg, obj_id; penalty_fn=alexs_penalty)\n\nPerform the teasar algorithm on the passed binary array.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.alexs_penalty",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.alexs_penalty",
    "category": "function",
    "text": "alexs_penalty( weights, dbf )\n\nReturns a version of the edge weights, modified by the DBF-based   teasar penalty:\n\nw = w * 5000 .* (1 - dbf/maxDBF).^(16)\n\nThe factor of 5000 shouldn\'t make a difference, but I\'ve   left it here\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.consolidate_paths-Tuple{Array{T,1} where T}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.consolidate_paths",
    "category": "method",
    "text": "consolidate_paths( path_list )\n\nExtracts the unique nodes and edges among all paths in the list\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.create_node_lookup-Union{Tuple{Array{T,2}}, Tuple{T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.create_node_lookup",
    "category": "method",
    "text": "create_node_lookup( points )\n\nAbstractly represent the points as a volume by means of   linear indexes into a sparse vector.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.distill!-Union{Tuple{T}, Tuple{Array{T,2},Array{T,1} where T,Array{T,1} where T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.distill!",
    "category": "method",
    "text": "Parameters:\n\npoint_array: a Nx3 array recording all the voxel coordinates inside the object.\npath_nodes: the indexes of nodeNet voxels in the point_array \npath_edges: the index pairs of nodeNet in the point_array\n\nIn the pathnodes, the nodes were encoded as the index of the pointarray. Since we do  not need pointarray anymore, which could be pretty big, we\'ll only reserve the nodeNet  coordinates and let the gc release the pointarray.  the same applys to path_edges \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.find_new_root_node_V1-Union{Tuple{Array{T,2}}, Tuple{T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.find_new_root_node_V1",
    "category": "method",
    "text": "find_new_root_node( points )\n\nExtracts the point with the lowest linear index from the   Array of passed points \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.find_new_root_node_id-Tuple{Array{Float32,1},DataStructures.IntSet}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.find_new_root_node_id",
    "category": "method",
    "text": "find_new_root_node( dbf )\n\nExtracts the point with the largest dbf from the    Array of passed points \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_connectivity_matrix-Tuple{Array{T,1} where T}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.get_connectivity_matrix",
    "category": "method",
    "text": "get_connectivity_matrix(edges::Vector)\n\nconstruct sparse connectivity matrix accordint to the edges \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_neuroglancer_precomputed-Tuple{RealNeuralNetworks.NodeNets.NodeNet}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.get_neuroglancer_precomputed",
    "category": "method",
    "text": "get binary buffer formatted as neuroglancer nodeNet.\n\nBinary format\n\nUInt32: number of vertex\nUInt32: number of edges\nArray{Float32,2}: Nx3 array, xyz coordinates of vertex\nArray{UInt32,2}: Mx2 arrray, node index pair of edges\n\nreference:  https://github.com/seung-lab/neuroglancer/wiki/Skeletons\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_num_segment_point-Tuple{RealNeuralNetworks.NodeNets.NodeNet}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.get_num_segment_point",
    "category": "method",
    "text": "get_segment_point_num(self::NodeNet)\n\nget number of branching points \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_num_segmentes-Tuple{RealNeuralNetworks.NodeNets.NodeNet}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.get_num_segmentes",
    "category": "method",
    "text": "assume that the graph is acyclic, no loop.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.get_sholl_number-Tuple{RealNeuralNetworks.NodeNets.NodeNet,AbstractFloat}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.get_sholl_number",
    "category": "method",
    "text": "get_sholl_number(self::NodeNet, radius::AbstractFloat)\n\nget the number of points which is in neurite and incounters with a sphere centered on root node \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.local_max_multiplicative_penalty-Tuple{Any,Any,Any}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.local_max_multiplicative_penalty",
    "category": "method",
    "text": "local_max_multiplicative_penalty( weights, dbf, G )\n\nReturns a version of the edge weights, modified by the DBF-based   penalty\n\nw = w * (1 - DBFdest/DBFstar)\n\nWhere DBFstar is the maximum DBF value of an outwards neighbor for   each node\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.make_neighbor_graph-Union{Tuple{Array{T,2}}, Tuple{T}, Tuple{Array{T,2},Any}, Tuple{Array{T,2},Any,Any}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.make_neighbor_graph",
    "category": "method",
    "text": "make_neighbor_graph( points )\n\nCreates a LightGraphs graph from a point cloud, such that   each point has a node (indexed in the same order as the rows),   and each node is connected to its neighbors (by 26-connectivity).\n\nAlso returns a sparse matrix of weights for these edges equal to the   euclidean distance between the indices of each point. These can   be weighted or modified easily upon return.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.nodes_within_radius-Union{Tuple{T}, Tuple{Array{T,1},Any,Any,Array{T,1} where T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.nodes_within_radius",
    "category": "method",
    "text": "nodes_within_radius( sub, sub2node, r, max_dims )\n\nIdentifies the node indices within the subscript\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.remove_path_from_rns!",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.remove_path_from_rns!",
    "category": "function",
    "text": "remove_path_from_rns!( reachableNodeList, path, points, sub2node, dbf, max_dims, scale_param, const_param)\n\nIdentifies the nodes to remove from the reachable nodes within the graph.   Probably the ugliest function here, and should be optimized later\n\nTO OPTIMIZE\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.save-Tuple{RealNeuralNetworks.NodeNets.NodeNet,Integer,AbstractDict}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.save",
    "category": "method",
    "text": "save(self::NodeNet, cellId::UInt32, d_json::Associative, d_bin::Associative)\n\nsave nodeNet in google cloud storage for neuroglancer visualization the format is the same with meshes\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.save_edges-Tuple{RealNeuralNetworks.NodeNets.NodeNet,String}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.save_edges",
    "category": "method",
    "text": "save binary file of point pairs used in neuroglancer python interface to visualize the nodeNet \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.translate_to_origin!-Tuple{Any}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.translate_to_origin!",
    "category": "method",
    "text": "translate_to_origin!( points )\n\nNormalize the point dimensions by subtracting the min   across each dimension. This step isn\'t extremely necessary,   but might be useful for compatibility with the MATLAB code.   record the offset and add it back after building the nodeNet\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{Array{T,2}}, Tuple{T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.compute_DBF",
    "category": "method",
    "text": "compute_DBF( pointCloud )\n\nReturns an array of DBF values for the point cloud. Currently creates   a binary image, and runs bwd2 on it, though ideally we\'d get rid of the   need for an explicit bin_im\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{T}, Tuple{Array{T,2},Array{T,1} where T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.compute_DBF",
    "category": "method",
    "text": "compute Distance from Boundary Field (DBF) based on point cloud and the boundary points\n\nWARN: this function do not work correctly!\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{T}, Tuple{Array{T,2},Union{Array{Bool,3}, BitArray}}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.compute_DBF",
    "category": "method",
    "text": "compute_DBF( bin_im )\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.compute_DBF-Union{Tuple{T}, Tuple{Array{T,3},T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.compute_DBF",
    "category": "method",
    "text": "use segmentation to get binary image to save memory usage\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.create_binary_image-Union{Tuple{Array{T,2}}, Tuple{T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.create_binary_image",
    "category": "method",
    "text": "create_binary_image( pointCloud )\n\nCreates a boolean volume where the non-segment indices   map to true, while the segment indices map to false.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.create_binary_image-Union{Tuple{Array{T,3}}, Tuple{T}} where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.create_binary_image",
    "category": "method",
    "text": "create_binary_image( seg, obj_id )\n\nCreates a boolean volume where the non-segment indices map to true, while the segment indices map to false \n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.distance_transform-Union{Tuple{AbstractArray{T,N}}, Tuple{N}, Tuple{T}, Tuple{AbstractArray{T,N},Array{Float32,1}}} where N where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.distance_transform",
    "category": "method",
    "text": "distance_transform( d::AbstractArray{T,N}, voxelSize::Vector{Float32}=ones(Float32, N) )\n\nReturns a euclidean distance transformation of the mask provided by d. The return   value will be a volume of the same size as d where the value at each index corresponds   to the distance between that location and the nearest location for which d > 0.\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.extract_dbf_values-Union{Tuple{N}, Tuple{Array{Float32,N},Any}} where N",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.extract_dbf_values",
    "category": "method",
    "text": "extract_dbf_values( dbf_image, pointCloud )\n\nTakes an array where rows indicate subscripts, and extracts the values   within a volume at those subscripts (in row order)\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.fill_f0!-Union{Tuple{N}, Tuple{T}, Tuple{Array{Float32,N},AbstractArray{T,N}}} where N where T",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.fill_f0!",
    "category": "method",
    "text": "Fills an n-dimensional volume with initial states for edt transformation, inf for non-feature voxels, and 0 for feature voxels\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.fv_isfurther-Tuple{Float32,Int64,Float32,Int64,Int64}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.fv_isfurther",
    "category": "method",
    "text": "Getting too tired to document these next few, but will be worth it if it works\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.remove_euclidean_distance_transform-Tuple{Float32,Float32,Float32,Int64,Int64,Int64}",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.remove_euclidean_distance_transform",
    "category": "method",
    "text": "remove_euclidean_distance_transform\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.row_voronoi_edt!-Union{Tuple{N}, Tuple{Array{Float32,N},Tuple,Array{Float32,1},Array{Int64,1},Float32}} where N",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.row_voronoi_edt!",
    "category": "method",
    "text": "Performs the edt over a specific row in the volume, following the first dimension\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#RealNeuralNetworks.NodeNets.DBFs.vol_voronoi_edt!-Union{Tuple{N}, Tuple{Array{Float32,N},Float32}} where N",
    "page": "NodeNets",
    "title": "RealNeuralNetworks.NodeNets.DBFs.vol_voronoi_edt!",
    "category": "method",
    "text": "Performs the edt transformation along the first dimension of the N-dimensional volume\n\n\n\n\n\n"
},

{
    "location": "lib/internals/NodeNets/#NodeNets-1",
    "page": "NodeNets",
    "title": "NodeNets",
    "category": "section",
    "text": "Modules = [\n    RealNeuralNetworks.NodeNets,\n    RealNeuralNetworks.NodeNets.DBFs\n]"
},

{
    "location": "lib/internals/PointArrays/#",
    "page": "PointArrays",
    "title": "PointArrays",
    "category": "page",
    "text": ""
},

{
    "location": "lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.add_offset!-Union{Tuple{T}, Tuple{Array{T,2},Tuple{T,T,T}}} where T",
    "page": "PointArrays",
    "title": "RealNeuralNetworks.NodeNets.PointArrays.add_offset!",
    "category": "method",
    "text": "add offset to points\n\n\n\n\n\n"
},

{
    "location": "lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.from_binary_image-Tuple{Array{Bool,3}}",
    "page": "PointArrays",
    "title": "RealNeuralNetworks.NodeNets.PointArrays.from_binary_image",
    "category": "method",
    "text": "parameter:     bin_im: binary array. object voxels are false, non-object voxels are true!\n\n\n\n\n\n"
},

{
    "location": "lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.from_seg-Union{Tuple{Array{T,3}}, Tuple{T}} where T",
    "page": "PointArrays",
    "title": "RealNeuralNetworks.NodeNets.PointArrays.from_seg",
    "category": "method",
    "text": "find points inside an object from a segmentation array. \n\n\n\n\n\n"
},

{
    "location": "lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.get_boundary_point_indexes-Union{Tuple{TSeg}, Tuple{T}, Tuple{Array{T,2},Array{TSeg,3}}, Tuple{Array{T,2},Array{TSeg,3},TSeg}} where TSeg where T",
    "page": "PointArrays",
    "title": "RealNeuralNetworks.NodeNets.PointArrays.get_boundary_point_indexes",
    "category": "method",
    "text": "find out the boundary voxels and represent them as indexes in the point array\n\n\n\n\n\n"
},

{
    "location": "lib/internals/PointArrays/#RealNeuralNetworks.NodeNets.PointArrays.merge-Union{Tuple{T}, Tuple{Array{T,2},Array{T,2}}} where T",
    "page": "PointArrays",
    "title": "RealNeuralNetworks.NodeNets.PointArrays.merge",
    "category": "method",
    "text": "merge(self::Array{T,3}, other::Array{T,3})\n\n\n\n\n\n"
},

{
    "location": "lib/internals/PointArrays/#PointArrays-1",
    "page": "PointArrays",
    "title": "PointArrays",
    "category": "section",
    "text": "Modules = [RealNeuralNetworks.NodeNets.PointArrays]"
},

{
    "location": "man/contributing/#",
    "page": "-",
    "title": "-",
    "category": "page",
    "text": "Please fork and create pull requests in github."
},

]}
