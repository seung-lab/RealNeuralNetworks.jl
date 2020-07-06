using Documenter
using RealNeuralNetworks
using RealNeuralNetworks.Neurons 
using RealNeuralNetworks.Neurons.Segments 
using RealNeuralNetworks.Neurons.Segments.Synapses
using RealNeuralNetworks.SWCs

# The DOCSARGS environment variable can be used to pass additional arguments to make.jl.
# This is useful on CI, if you need to change the behavior of the build slightly but you
# can not change the .travis.yml or make.jl scripts any more (e.g. for a tag build).
if haskey(ENV, "DOCSARGS")
    for arg in split(ENV["DOCSARGS"])
        (arg in ARGS) || push!(ARGS, arg)
    end
end

makedocs(
	modules=[RealNeuralNetworks, Neurons, Segments, Synapses, SWCs],
	sitename="RealNeuralNetworks.jl",
	authors="Jingpeng Wu",
    format = Documenter.HTML(
        # Use clean URLs, unless built as a "local" build
        prettyurls = !("local" in ARGS),
        canonical = "https://seung-lab.github.io/RealNeuralNetworks.jl/latest/",
        assets = ["assets/favicon.ico"],
        analytics = "UA-136089579-2",
        highlights = ["yaml"],
    ),
    linkcheck = !("skiplinks" in ARGS),
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
)

deploydocs(
	repo="github.com/seung-lab/RealNeuralNetworks.jl",
	target="build",
	deps=nothing,
	make=nothing
)
