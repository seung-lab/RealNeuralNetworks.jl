using Documenter, RealNeuralNetworks

makedocs(
	modules=[RealNeuralNetworks],
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
	julia="1.0.0",
	deps=nothing,
	make=nothing
)
