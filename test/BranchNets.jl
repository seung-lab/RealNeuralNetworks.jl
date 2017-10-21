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
    println("create fake cylinder segmentation...")
    @time seg = FakeSegmentations.broken_cylinder()
    println("skeletonization to build a BranchNet ...")
    @time branchNet = BranchNet(seg)
    println("transform to SWC structure ...")
    @time swc = SWCs.SWC( branchNet )
    SWCs.save(swc, "/tmp/cylinder.swc")

    println("create fake ring segmentation ...")
    seg = FakeSegmentations.broken_ring()
    branchNet = BranchNet(seg)
    swc = SWCs.SWC( branchNet )
    SWCs.save(swc, "/tmp/ring.swc")

    println("get terminal branch index list...")
    @time terminalBranchIndexList = BranchNets.get_terminal_branch_index_list( branchNet )
    @test !isempty( terminalBranchIndexList )
    println("get terminal node list ...")
    @time terminalNodeList = BranchNets.get_terminal_node_list( branchNet )
    @test !isempty( terminalNodeList )
end 
