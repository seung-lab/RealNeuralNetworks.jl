#!/usr/bin/env julia

using RealNeuralNetworks
using RealNeuralNetworks.SWCs
using RealNeuralNetworks.BranchNets 

include("Common.jl"); using Common

const CELL_ID = 77390

function main()
    args = parse_commandline()
    swc = SWCs.load_gzip_swc("swcgz/$(CELL_ID).swc.gz")
    branchNet = BranchNet( swc )
    branchNet = BranchNets.remove_subtree_in_soma(branchNet)
    branchNet = BranchNets.remove_hair( branchNet )
    BranchNets.save(branchNet, "$(CELL_ID).swc")
end

main()
