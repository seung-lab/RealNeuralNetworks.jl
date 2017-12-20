#!/usr/bin/env julia
using Plots
pyplot()

using RealNeuralNetworks.Neurons
import RealNeuralNetworks.Neurons.Segments


neuron = Neurons.load_swc("/tmp/76391.swc")
segmentList = Neurons.get_segment_list(neuron)
global x,y,z

for segment in segmentList
    nodeList = Segments.get_node_list(segment)
    #@show length(nodeList)
    #@show nodeList
    x = map(n->n[1], nodeList)
    y = map(n->n[2], nodeList)
    z = map(n->n[3], nodeList)
    plot!(x,y,z) #color=rand(Colors.RGB))
end

root = Neurons.get_root_node( neuron )
display(plot!([root[1]],[root[2]],[root[3]], m=(8,:auto), leg=true))

println("start sleeping ...")
sleep(3000)
