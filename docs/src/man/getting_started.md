
```Julia
using RealNeuralNetworks.NodeNets
using RealNeuralNetworks.Neurons
using RealNeuralNetworks.SWCs

# skeletonization
nodeNet = NodeNet(seg::Array{UInt32,3}; obj_id = convert(UInt32,77605))
neuron = Neuron( nodeNet )
swc = SWC(neuron)
SWCs.save(swc, tempname()*".swc")
```
