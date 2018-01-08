__precompile__();
module RealNeuralNetworks

include("FakeSegmentations.jl")
include("SWCs.jl");  
include("NodeNets.jl"); 
include("Manifests.jl"); 
include("Neurons.jl")
include("NeuralNets.jl")

using .SWCs
using .NodeNets 
using .Manifests
using .Neurons
using .NeuralNets

end # end of module
