using Test
using RealNeuralNetworks.Utils.FakeSegmentations 
using RealNeuralNetworks.Neurons
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs
using RealNeuralNetworks.Utils.SynapseTables 

using CSV
using SparseArrays 
using LinearAlgebra 
using DataFrames 

const NEURON_ID = 77625
const ASSET_DIR = joinpath(@__DIR__, "../asset")
const SWC_BIN_PATH = joinpath(ASSET_DIR, "$(NEURON_ID).swc.bin") 
const ARBOR_DENSITY_MAP_VOXEL_SIZE = (2000,2000,2000)

@testset "test data structure transformation..." begin 
    swc = SWCs.load(joinpath(ASSET_DIR, "example.swc"))
    nodeNet = NodeNet(swc)
    neuron = Neuron(nodeNet)
    swc = SWC(neuron)
    nodeNet = NodeNet(neuron)
    @test length(swc) == length(nodeNet) == Neurons.get_num_nodes(neuron)
end 

@testset "test Neurons" begin
    println("\nload swc of a real neuron...")
    swc = SWCs.load_swc_bin( SWC_BIN_PATH )
    neuron = Neuron( swc )

    println("\nreset root node...")
    @time Neurons.reset_root(neuron, (258150.0f0, 135510.0f0, 778175.0f0))

    #neuron = Neurons.resample(neuron, Float32(40))
    println("\nget node list ...")
    @time nodeList = Neurons.get_node_list(neuron)
    @test !isempty(nodeList) 
    println("\nget edge list ...")
    @time edgeList = Neurons.get_edge_list(neuron)
    @test !isempty(edgeList)
    println("\nget segment order list...")
    @time segmentOrderList = Neurons.get_segment_order_list( neuron )
    @test !isempty(segmentOrderList) 

    println("\nclean up the neuron ...")
    println("original number of segments: $(Neurons.get_num_segments(neuron))")
    println("remove subtree in soma ...")
    @time neuron = Neurons.remove_subtree_in_soma(neuron)
    println("after remove subtree in soma: $(Neurons.get_num_segments(neuron))")

    println("\nremove hair segments ...")
    @time neuron = Neurons.remove_hair(neuron)
    println("after remove hair: $(Neurons.get_num_segments(neuron))")
 
    println("\nremove terminal blobs ...")
    @time neuron = Neurons.remove_terminal_blobs(neuron)
    println("after remove terminal blobs: $(Neurons.get_num_segments(neuron))") 

    println("\nget principle direction based on the furthest terminal node pair...")
    @time vec, maxDistance = Neurons.get_furthest_terminal_node_pair_direction(neuron)
    @test maxDistance > 0
  
    println("\nremove redundent nodes...")
    @time neuron = Neurons.remove_redundent_nodes(neuron)
    println("after remove redundent nodes: $(Neurons.get_num_segments(neuron))")

   
    println("\nresampling ...")
    nodeNum1 = Neurons.get_node_num(neuron)
    @time neuron = Neurons.resample(neuron, Float32(400))
    nodeNum2 = Neurons.get_node_num(neuron)
    println("resampling of 400 nm reduced node number form ", nodeNum1, " to ", nodeNum2)
    @test nodeNum2 < nodeNum1
 
    println("\nsmooth ...")
    nodeNum1 = Neurons.get_node_num(neuron) 
    @time neuron = Neurons.smooth(neuron)
    nodeNum2 = Neurons.get_node_num(neuron)
    println("smoothing the neuron should not change node number: ", nodeNum1, " => ", nodeNum2)
    @test nodeNum1 == nodeNum2

    
    @test Neurons.get_num_segments(neuron) > 0
    #println("get fractal dimension ...")
    #@time fractalDimension, _,_ = Neurons.get_fractal_dimension( neuron )
    #@show fractalDimension 

    println("\nget surface area which is frustum based")
    @test Neurons.get_surface_area(neuron) > 0
    
    println("\nget frustum based volume")
    @test Neurons.get_volume(neuron) > 0

    println("\nget typical radius ...")
    @test Neurons.get_typical_radius( neuron ) > 0

    println("\nget asymmetry ...")
    @test Neurons.get_asymmetry( neuron ) > 0
 
    println("\nget mass center ...")
    @show Neurons.get_mass_center( neuron )

    println("\nget branching angle ...")
    @test Neurons.get_branching_angle( neuron, 5 ) > 0

    println("\nget path to root length ...")
    @test Neurons.get_path_to_soma_length( neuron, 5; nodeId=4 ) > 0

    println("\nsholl analysis ...")
    @time shollNumList = Neurons.get_sholl_number_list(neuron, 10000 )
    @test !isempty(shollNumList)

    println("\nget segment path length list ...")
    @time segmentPathLengthList = Neurons.get_segment_path_length_list( neuron )
    @show length( segmentPathLengthList )
    @test length( segmentPathLengthList ) == Neurons.get_num_segments(neuron)

    println("\nget terminal segment index list...")
    @time terminalSegmentIdList = Neurons.get_terminal_segment_id_list( neuron )
    @test !isempty( terminalSegmentIdList )
    
    println("\nget terminal node list ...")
    @time terminalNodeList = Neurons.get_terminal_node_list( neuron )
    @test !isempty( terminalNodeList )

    println("\ntest split a tree...")
    @time tree1, tree2 = split(neuron, 4; nodeIdInSegment = 5)
    @test !isempty(tree1)
    @test !isempty(tree2)
end 


@testset "test synapse attaching functions" begin
    println("\nload swc of a real neuron...")
    @time swc = SWCs.load_swc_bin( SWC_BIN_PATH )
    @time neuron = Neuron( swc )

    println("\ndownsample node number...")
    @time Neurons.downsample_nodes(neuron)

    println("\nattaching presynapses in DataFrame...")
    # the default reading of CSV is immutable and the dropmissing will fail!
    # https://github.com/JuliaData/DataFrames.jl/issues/1393#issuecomment-500019773
    preSynapses = CSV.read( joinpath(ASSET_DIR, "$(NEURON_ID).pre.synapses.csv");  copycols=true)
    @time preSynapses = SynapseTables.preprocess(preSynapses, (5,5,45))
    println("get ", DataFrames.nrow(preSynapses), " synapses.")
    @time Neurons.attach_pre_synapses!(neuron, preSynapses)
    numPreSynapses = Neurons.get_num_pre_synapses(neuron)
    println("have ", numPreSynapses, " synapses after attachment.")
    @test numPreSynapses <= DataFrames.nrow(preSynapses) 
    @test numPreSynapses > 2 
    
    println("\nattaching postsynapses in DataFrame...")
    postSynapses = CSV.read( joinpath(ASSET_DIR, "$(NEURON_ID).post.synapses.csv"); copycols=true) 
    @time postSynapses = SynapseTables.preprocess(postSynapses, (5,5,45))
    println("get ", DataFrames.nrow(postSynapses), " synapses.")
    @time Neurons.attach_post_synapses!(neuron, postSynapses)
    numPostSynapses = Neurons.get_num_post_synapses(neuron)
    println("have ", numPostSynapses, " synapses after attachment.")
    @test numPostSynapses <= DataFrames.nrow(postSynapses) 
    @test numPostSynapses > 2 

    println("\nattaching presynapse list...")
    preSynapseList = Neurons.get_all_pre_synapse_list(neuron)
    println("get ", length(preSynapseList), " synapses.")
    neuronx = Neuron(swc)
    @time Neurons.attach_pre_synapses!(neuronx, preSynapseList)
    numPreSynapses = Neurons.get_num_pre_synapses(neuronx)
    println("have ", numPreSynapses, " synapses after attachment.")
    @test numPreSynapses > 2

    println("\nattaching postsynapse list...")
    postSynapseList = Neurons.get_all_post_synapse_list(neuron)
    println("get ", length(postSynapseList), " synapses.")
    @time Neurons.attach_post_synapses!(neuronx, postSynapseList)
    numPostSynapses = Neurons.get_num_post_synapses(neuronx)
    println("have ", numPostSynapses, " synapses after attachment.")
    @test numPostSynapses > 2

    println("\nresample neuron and double check the synapses.")
    neuron2 = Neurons.resample(neuron, Float32(1000))
    @time numPostSynapses2 = Neurons.get_num_post_synapses(neuron2)
    @time numPreSynapses2 = Neurons.get_num_pre_synapses(neuron2)
    @test numPostSynapses >= numPostSynapses2
    @test numPreSynapses >= numPreSynapses2 

    println("\npostprocessing neuron and make sure that the synapses are still there.")
    neuron2 = Neurons.postprocessing(neuron)
    @time numPostSynapses2 = Neurons.get_num_post_synapses(neuron2)
    @time numPreSynapses2 = Neurons.get_num_pre_synapses(neuron2)
    @test numPostSynapses >= numPostSynapses2
    @test numPreSynapses >= numPreSynapses2

    
    println("\nadjust segment class...")
    @time Neurons.adjust_segment_class!(neuron)
    
    preSynapseList = Neurons.get_all_pre_synapse_list(neuron)
    postSynapseList = Neurons.get_all_post_synapse_list(neuron)
    @test !isempty(preSynapseList) 
    @test !isempty(postSynapseList)

    println("\nget synapse to soma path length list...")
    @time preSynapseToSomaPathLengthList  = Neurons.get_pre_synapse_to_soma_path_length_list( neuron )
    @time postSynapseToSomaPathLengthList = Neurons.get_post_synapse_to_soma_path_length_list( neuron )
    @test !isempty(preSynapseToSomaPathLengthList)
    @test !isempty(postSynapseToSomaPathLengthList)
end 

@testset "test Neuron IO and resampling " begin 
    println("\nload swc of a real neuron...")
    @time swc = SWCs.load_swc_bin( SWC_BIN_PATH )
    neuron = Neuron( swc )

    println("\ndownsample node number...")
    Neurons.downsample_nodes(neuron)

    println("\ncompute arbor density map...")
    @time arborDensityMap = Neurons.get_arbor_density_map(neuron, 
                                                ARBOR_DENSITY_MAP_VOXEL_SIZE, 8.0)
    #@test norm(arborDensityMap[:]) ≈ Neurons.get_total_path_length(neuron)
    #@test norm(arborDensityMap[:]) ≈ 1.0
    neuron1 = neuron
    densityMap1 = Neurons.translate_soma_to_coordinate_origin(neuron1, arborDensityMap, 
                                                             ARBOR_DENSITY_MAP_VOXEL_SIZE)
    @show axes(densityMap1)
    println("compute arbor density map distance...")
    @time d = Neurons.get_arbor_density_map_distance(densityMap1, densityMap1)
    @test d == 0.0

    fileName = joinpath(dirname(SWC_BIN_PATH), "$(NEURON_ID).swc.bin")
    neuron2 = Neuron( SWCs.load_swc_bin(fileName) )
    densityMap2 = Neurons.get_arbor_density_map(neuron, ARBOR_DENSITY_MAP_VOXEL_SIZE, 8.0)
    #@test norm(densityMap2[:]) ≈ 1.0                                               
    densityMap2 = Neurons.translate_soma_to_coordinate_origin(neuron2, densityMap2, 
                                                              ARBOR_DENSITY_MAP_VOXEL_SIZE)
    @show axes(densityMap2) 
    println("compute arbor density map distance...")                          
    @time d = Neurons.get_arbor_density_map_distance(densityMap1, densityMap2)
    println("arbor density map distance: $d")
    #@test d > 0.0 && d < 2.0

    Neurons.save(neuron, "/tmp/neuron.swc")
    #Neurons.save_swc(neuron2, "/tmp/neuron2.swc")
    rm("/tmp/neuron.swc")
    #rm("/tmp/neuron2.swc")
end 

@testset "test fake segmentation skeletonization" begin 
    println("create fake cylinder segmentation...")
    @time seg = FakeSegmentations.broken_cylinder()
    println("skeletonization to build a Neuron ...")
    @time neuron = Neuron(seg)
    @test !isempty(Neurons.get_node_list(neuron))
    
    println("create fake ring segmentation ...")
    seg = FakeSegmentations.broken_ring()
    neuron = Neuron(seg)
    @test !isempty(Neurons.get_node_list(neuron))
    
    println("transform to SWC structure ...")
    @time swc = SWCs.SWC( neuron )
    tempFile = tempname() * ".swc"
    SWCs.save(swc, tempFile)
    rm(tempFile)
end 


