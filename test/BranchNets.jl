using Base.Test
using RealNeuralNetworks.FakeSegmentations 
using RealNeuralNetworks.BranchNets
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.SWCs

@testset "test BranchNet IO" begin 
    branchNet = BranchNets.load_swc( joinpath(dirname(@__FILE__), "../assert/example.swc" ))
    BranchNets.save(branchNet, "/tmp/branchNet.swc")
    rm("/tmp/branchNet.swc")
end 

@testset "test BranchNets" begin
    #println("create fake cylinder segmentation...")
    #@time seg = FakeSegmentations.broken_cylinder()
    #println("skeletonization to build a BranchNet ...")
    #@time branchNet = BranchNet(seg)
    println("create fake ring segmentation ...")
    seg = FakeSegmentations.broken_ring()
    branchNet = BranchNet(seg)
    println("transform to SWC structure ...")
    @time swc = SWCs.SWC( branchNet )
    tempFile = tempname() * ".swc"
    SWCs.save(swc, tempFile)
    rm(tempFile)

    println("load swc of a real neuron...")
    @time swc = SWCs.load_gzip_swc("../assert/76869.swc.gz")
    println("save as gzip compressed file ...")
    @time SWCs.save_gzip_swc(swc, "/tmp/76869.swc.gz")
    rm("/tmp/76869.swc.gz")

    branchNet = BranchNet( swc )
    println("get node list ...")
    @time nodeList = BranchNets.get_node_list(branchNet)
    println("get edge list ...")
    @time edgeList = BranchNets.get_edge_list(branchNet)

    println("remove subtree in soma...")
    @time newBranchNet = BranchNets.remove_subtree_in_soma(branchNet)

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
