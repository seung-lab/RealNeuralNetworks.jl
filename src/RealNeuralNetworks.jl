module RealNeuralNetworks

include("Utils/include.jl"); 
include("NodeNets.jl"); 
include("Manifests.jl"); 
include("Neurons.jl")
include("NeuralNets.jl")
include("NBLASTs.jl")

using .Manifests
using .Neurons
using .NeuralNets
using .NBLASTs

end # end of module
