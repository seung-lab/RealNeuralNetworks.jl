using Base.Test
using RealNeuralNetworks.FakeSegmentations 
using RealNeuralNetworks.BranchNets
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs

@testset "test BranchNet IO and resampling " begin 
    println("load swc of a real neuron...")
    @time swc = SWCs.load_swc_bin("../assert/78058.swc.bin")
    branchNet = BranchNet( swc )
    
    BranchNets.save(branchNet, "/tmp/branchNet.swc")
    branchNet2 = BranchNets.resample(branchNet, Float32(10))
    BranchNets.save_swc(branchNet2, "/tmp/branchNet2.swc")
    #rm("/tmp/branchNet.swc")
end 

@testset "test BranchNets" begin
    println("load swc of a real neuron...")
    @time swc = SWCs.load_swc_bin("../assert/78058.swc.bin")

    branchNet = BranchNet( swc )
    println("get node list ...")
    @time nodeList = BranchNets.get_node_list(branchNet)
    println("get edge list ...")
    @time edgeList = BranchNets.get_edge_list(branchNet)
    println("get branch order list...")
    @time branchOrderList = BranchNets.get_branch_order_list( branchNet )

    println("clean up the neuron ...")
    branchNet = BranchNets.remove_subtree_in_soma(branchNet)
    branchNet = BranchNets.remove_hair(branchNet)
    branchNet = BranchNets.remove_subtree_in_soma(branchNet)
    branchNet = BranchNets.remove_terminal_blobs(branchNet)
    branchNet = BranchNets.remove_redundent_nodes(branchNet)

    println("get fractal dimension ...")
    @time fractalDimension, _,_ = BranchNets.get_fractal_dimension( branchNet )
    @show fractalDimension 

    println("get typical radius ...")
    @show BranchNets.get_typical_radius( branchNet )

    println("get asymmetry ...")
    @show BranchNets.get_asymmetry( branchNet )
 
    println("get mass center ...")
    @show BranchNets.get_mass_center( branchNet )

    println("get branching angle ...")
    @time angle = BranchNets.get_branching_angle( branchNet, 5 )

    println("get path to root length ...")
    @time path2RootLength = BranchNets.get_path_to_root_length( branchNet, 5 )

    println("sholl analysis ...")
    @time shollNumList = BranchNets.get_sholl_number_list(branchNet, 10000 )

    println("get branch path length list ...")
    @time branchPathLengthList = BranchNets.get_branch_path_length_list( branchNet )
    @show branchPathLengthList
    @show length( branchPathLengthList )
    @show BranchNets.get_num_branches( branchNet ) 
    @test length( branchPathLengthList ) == BranchNets.get_num_branches(branchNet)

    println("get terminal branch index list...")
    @time terminalBranchIndexList = BranchNets.get_terminal_branch_index_list( branchNet )
    @test !isempty( terminalBranchIndexList )
    println("get terminal node list ...")
    @time terminalNodeList = BranchNets.get_terminal_node_list( branchNet )
    @test !isempty( terminalNodeList )

    println("test split a tree...")
    @time tree1, tree2 = split(branchNet, 2; nodeIndexInBranch = 5)
    @test !isempty(tree1)
    @test !isempty(tree2)
end 

@testset "test fake segmentation skeletonization" begin 
    println("create fake cylinder segmentation...")
    @time seg = FakeSegmentations.broken_cylinder()
    println("skeletonization to build a BranchNet ...")
    @time branchNet = BranchNet(seg)
    
    println("create fake ring segmentation ...")
    seg = FakeSegmentations.broken_ring()
    branchNet = BranchNet(seg)
    println("transform to SWC structure ...")
    @time swc = SWCs.SWC( branchNet )
    tempFile = tempname() * ".swc"
    SWCs.save(swc, tempFile)
    rm(tempFile)
end 


