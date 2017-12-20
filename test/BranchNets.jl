using Base.Test
using RealNeuralNetworks.FakeSegmentations 
using RealNeuralNetworks.BranchNets
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs

@testset "test BranchNet IO and resampling " begin 
    println("load swc of a real neuron...")
    @time swc = SWCs.load_swc_bin("../assert/78058.swc.bin")
    neuron = BranchNet( swc )
    
    arborDensityMap = BranchNets.get_arbor_density_map(neuron, (360,360,360), 8.0)
    @test norm(arborDensityMap[:]) ≈ BranchNets.get_total_path_length(neuron)
    
    BranchNets.save(neuron, "/tmp/neuron.swc")
    neuron2 = BranchNets.resample(neuron, Float32(40))
    BranchNets.save_swc(neuron2, "/tmp/neuron2.swc")
    rm("/tmp/neuron.swc")
    rm("/tmp/neuron2.swc")
end 

@testset "test BranchNets" begin
    println("load swc of a real neuron...")
    @time swc = SWCs.load_swc_bin("../assert/78058.swc.bin")

    neuron = BranchNet( swc )
    println("get node list ...")
    @time nodeList = BranchNets.get_node_list(neuron)
    println("get edge list ...")
    @time edgeList = BranchNets.get_edge_list(neuron)
    println("get branch order list...")
    @time branchOrderList = BranchNets.get_branch_order_list( neuron )

    println("clean up the neuron ...")
    neuron = BranchNets.remove_subtree_in_soma(neuron)
    neuron = BranchNets.remove_hair(neuron)
    neuron = BranchNets.remove_subtree_in_soma(neuron)
    neuron = BranchNets.remove_terminal_blobs(neuron)
    neuron = BranchNets.remove_redundent_nodes(neuron)

    println("get fractal dimension ...")
    @time fractalDimension, _,_ = BranchNets.get_fractal_dimension( neuron )
    @show fractalDimension 

    println("get typical radius ...")
    @show BranchNets.get_typical_radius( neuron )

    println("get asymmetry ...")
    @show BranchNets.get_asymmetry( neuron )
 
    println("get mass center ...")
    @show BranchNets.get_mass_center( neuron )

    println("get branching angle ...")
    @time angle = BranchNets.get_branching_angle( neuron, 5 )

    println("get path to root length ...")
    @time path2RootLength = BranchNets.get_path_to_root_length( neuron, 5 )

    println("sholl analysis ...")
    @time shollNumList = BranchNets.get_sholl_number_list(neuron, 10000 )

    println("get branch path length list ...")
    @time branchPathLengthList = BranchNets.get_branch_path_length_list( neuron )
    @show branchPathLengthList
    @show length( branchPathLengthList )
    @show BranchNets.get_num_branches( neuron ) 
    @test length( branchPathLengthList ) == BranchNets.get_num_branches(neuron)

    println("get terminal branch index list...")
    @time terminalBranchIndexList = BranchNets.get_terminal_branch_index_list( neuron )
    @test !isempty( terminalBranchIndexList )
    println("get terminal node list ...")
    @time terminalNodeList = BranchNets.get_terminal_node_list( neuron )
    @test !isempty( terminalNodeList )

    println("test split a tree...")
    @time tree1, tree2 = split(neuron, 2; nodeIndexInBranch = 5)
    @test !isempty(tree1)
    @test !isempty(tree2)
end 

@testset "test fake segmentation skeletonization" begin 
    println("create fake cylinder segmentation...")
    @time seg = FakeSegmentations.broken_cylinder()
    println("skeletonization to build a BranchNet ...")
    @time neuron = BranchNet(seg)
    
    println("create fake ring segmentation ...")
    seg = FakeSegmentations.broken_ring()
    neuron = BranchNet(seg)
    println("transform to SWC structure ...")
    @time swc = SWCs.SWC( neuron )
    tempFile = tempname() * ".swc"
    SWCs.save(swc, tempFile)
    rm(tempFile)
end 


