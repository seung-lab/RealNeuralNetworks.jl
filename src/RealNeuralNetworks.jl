__precompile__();
module RealNeuralNetworks

include("FakeSegmentations.jl")
include("SWCs.jl");  
include("NodeNets.jl"); 
include("Manifests.jl"); 
include("Neurons.jl")

using .SWCs
using .NodeNets 
using .Manifests
using .Neurons

end # end of module
