using RealNeuralNetworks
using RealNeuralNetworks.BranchNets
using RealNeuralNetworks.BranchNets.Branches
using Base.Test


@testset "test Branches" begin
    # construct a branch
    branchNet = BranchNets.load_swc_bin("../assert/78058.swc.bin")
    # get a random branch
    branch = branchNet[5]
    
    println("get tortuosity...")
    @show Branches.get_tortuosity( branch )
    println("get tail head radius ratio ...")
    @show Branches.get_tail_head_radius_ratio( branch )
end 
