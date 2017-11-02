#!/usr/bin/env julia
using Plots
pyplot()

using RealNeuralNetworks.BranchNets
import RealNeuralNetworks.BranchNets.Branches


branchNet = BranchNets.load_swc("/tmp/76391.swc")
branchList = BranchNets.get_branch_list(branchNet)
global x,y,z

for branch in branchList
    nodeList = Branches.get_node_list(branch)
    #@show length(nodeList)
    #@show nodeList
    x = map(n->n[1], nodeList)
    y = map(n->n[2], nodeList)
    z = map(n->n[3], nodeList)
    plot!(x,y,z) #color=rand(Colors.RGB))
end

root = BranchNets.get_root_node( branchNet )
display(plot!([root[1]],[root[2]],[root[3]], m=(8,:auto), leg=true))

println("start sleeping ...")
sleep(3000)
