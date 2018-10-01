using Documenter, RealNeuralNetworks
using RealNeuralNetworks.Neurons 
using RealNeuralNetworks.Neurons.Segments 
using RealNeuralNetworks.Neurons.Segments.Synapses
using RealNeuralNetworks.SWCs

makedocs(
	modules=[RealNeuralNetworks, Neurons, Segments, Synapses, SWCs],
	sitename="RealNeuralNetworks.jl",
	authors="Jingpeng Wu",
	format=:html,
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Getting Started" => "man/getting_started.md",
            "man/examples.md",
            "man/references.md",
        ],
        "Library" => Any[
            "Public" => "lib/public.md",
            hide("Internals" => "lib/internals.md", Any[
                "lib/internals/Manifests.md",
                "lib/internals/NodeNets.md",
                "lib/internals/PointArrays.md",
            ])
        ],
        "man/contributing.md",
    ],
    # use clean URLs, unless built as a "local" build 
    html_prettyurls = !("local" in ARGS),
    html_canonical = "https://seung-lab.github.io/RealNeuralNetworks.jl/latest",
)

deploydocs(
	repo="github.com/seung-lab/RealNeuralNetworks.jl",
	target="build",
	julia="0.7",
	deps=nothing,
	make=nothing
)
