using Documenter, RealNeuralNetworks
using RealNeuralNetworks.Neurons 
using RealNeuralNetworks.Neurons.Segments 
using RealNeuralNetworks.Neurons.Segments.Synapses
using RealNeuralNetworks.SWCs

makedocs(
	modules=[RealNeuralNetworks, Neurons, Segments, Synapses, SWCs],
	sitename="RealNeuralNetworks.jl",
	authors="Jingpeng Wu",
	format=:html
)

deploydocs(
	repo="github.com/seung-lab/RealNeuralNetworks.jl",
	target="build",
	julia="0.7",
	deps=nothing,
	make=nothing
)
