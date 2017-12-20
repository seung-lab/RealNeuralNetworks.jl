using Base.Test
using RealNeuralNetworks.FakeSegmentations 
using RealNeuralNetworks.Neurons
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs

const SWC_BIN_PATH = joinpath(@__DIR__, "../assert/78058.swc.bin") 

@testset "test Neuron IO and resampling " begin 
    println("load swc of a real neuron...")
    @time swc = SWCs.load_swc_bin( SWC_BIN_PATH )
    neuron = Neuron( swc )
    
    arborDensityMap = Neurons.get_arbor_density_map(neuron, (360,360,360), 8.0)
    @test norm(arborDensityMap[:]) ≈ Neurons.get_total_path_length(neuron)
    
    Neurons.save(neuron, "/tmp/neuron.swc")
    neuron2 = Neurons.resample(neuron, Float32(40))
    #Neurons.save_swc(neuron2, "/tmp/neuron2.swc")
    rm("/tmp/neuron.swc")
    #rm("/tmp/neuron2.swc")
end 

@testset "test Neurons" begin
    println("load swc of a real neuron...")
    @time swc = SWCs.load_swc_bin( SWC_BIN_PATH )

    neuron = Neuron( swc )
    #neuron = Neurons.resample(neuron, Float32(40))
    println("get node list ...")
    @time nodeList = Neurons.get_node_list(neuron)
    println("get edge list ...")
    @time edgeList = Neurons.get_edge_list(neuron)
    println("get segment order list...")
    @time segmentOrderList = Neurons.get_segment_order_list( neuron )

    println("clean up the neuron ...")
    neuron = Neurons.remove_subtree_in_soma(neuron)
    neuron = Neurons.remove_hair(neuron)
    neuron = Neurons.remove_subtree_in_soma(neuron)
    neuron = Neurons.remove_terminal_blobs(neuron)
    neuron = Neurons.remove_redundent_nodes(neuron)

    println("get fractal dimension ...")
    @time fractalDimension, _,_ = Neurons.get_fractal_dimension( neuron )
    @show fractalDimension 

    println("get typical radius ...")
    @show Neurons.get_typical_radius( neuron )

    println("get asymmetry ...")
    @show Neurons.get_asymmetry( neuron )
 
    println("get mass center ...")
    @show Neurons.get_mass_center( neuron )

    println("get segmenting angle ...")
    @time angle = Neurons.get_segmenting_angle( neuron, 5 )

    println("get path to root length ...")
    @time path2RootLength = Neurons.get_path_to_root_length( neuron, 5 )

    println("sholl analysis ...")
    @time shollNumList = Neurons.get_sholl_number_list(neuron, 10000 )

    println("get segment path length list ...")
    @time segmentPathLengthList = Neurons.get_segment_path_length_list( neuron )
    @show segmentPathLengthList
    @show length( segmentPathLengthList )
    @show Neurons.get_num_segmentes( neuron ) 
    @test length( segmentPathLengthList ) == Neurons.get_num_segmentes(neuron)

    println("get terminal segment index list...")
    @time terminalBranchIndexList = Neurons.get_terminal_segment_index_list( neuron )
    @test !isempty( terminalBranchIndexList )
    println("get terminal node list ...")
    @time terminalNodeList = Neurons.get_terminal_node_list( neuron )
    @test !isempty( terminalNodeList )

    println("test split a tree...")
    @time tree1, tree2 = split(neuron, 4; nodeIndexInBranch = 5)
    @test !isempty(tree1)
    @test !isempty(tree2)
end 

@testset "test fake segmentation skeletonization" begin 
    println("create fake cylinder segmentation...")
    @time seg = FakeSegmentations.broken_cylinder()
    println("skeletonization to build a Neuron ...")
    @time neuron = Neuron(seg)
    
    println("create fake ring segmentation ...")
    seg = FakeSegmentations.broken_ring()
    neuron = Neuron(seg)
    println("transform to SWC structure ...")
    @time swc = SWCs.SWC( neuron )
    tempFile = tempname() * ".swc"
    SWCs.save(swc, tempFile)
    rm(tempFile)
end 


