using RealNeuralNetworks
using RealNeuralNetworks.BranchNets
using RealNeuralNetworks.BranchNets.Branches
using Base.Test


@testset "test Branches" begin
    # construct a branch
    branchNet = BranchNets.load_gzip_swc("../assert/76869.swc.gz")
    # get a random branch
    branch = branchNet[5]
    
    println("get tail head radius ratio ...")
    @time Branches.get_tail_head_radius_ratio( branch )
end 
