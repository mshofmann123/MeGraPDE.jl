push!(LOAD_PATH, "../src/")

using Pkg
Pkg.instantiate()

using Documenter
using MeGraPDE

index_page = "index.md"

gettingstarted_page="Getting Started"=>[
    "gettingstarted/metricgraphs.md",
    "gettingstarted/spectrum.md",
    "gettingstarted/heat.md",
    "gettingstarted/theory.md",
]

metricgraphs_page =
    "Metric Graphs" => [
        "metricgraphs/basicusage.md",
        "metricgraphs/constructors.md",
        "metricgraphs/extendedgraph.md",
    ]

finite_elements_page="Finite Elements"=>["finite_elements/discretization.md"]
plots_page = "plots.md"
spectra_page = "Spectra"=>["spectra/equilateral.md", "spectra/non_equilateral.md"]
spectral_galerkin_page = "spectral_galerkin.md"
testproblems_page = "testproblems.md"

makedocs(;
    sitename="MeGraPDE",
    pages=[
        index_page,
        gettingstarted_page,
        metricgraphs_page,
        finite_elements_page,
        plots_page,
        spectra_page,
        spectral_galerkin_page,
        testproblems_page,
    ],
)
