push!(LOAD_PATH,"../src/")
using Documenter, MeGraPDE, DocStringExtensions
using Graphs, Plots
include("../src/finite_elements/multigrid_matrix_free.jl")

SUBSECTION_GettingStarted = ["gettingstarted/metricgraphs.md", "gettingstarted/spectrum.md", "gettingstarted/heat.md"]
SUBSECTION_MetricGraphs = ["metricgraphs/basicusage.md", "metricgraphs/extendedgraph.md", "metricgraphs/constructors.md"]
SUBSECTION_FiniteElements = ["finite_elements/discretization.md", "finite_elements/multigrid.md"]
SUBSECTION_Spectra = ["spectra/equilateral.md", "spectra/non_equilateral.md"]
ENV["PLOTS_TEST"] = "true" 
ENV["GKSwstype"] = "100"
makedocs(
    sitename="MeGraPDE",
    pages = [
        "index.md",
        "Getting Started" => SUBSECTION_GettingStarted,
        "Metric Graphs" => SUBSECTION_MetricGraphs,
        "testproblems.md",
        "Finite Elements" => SUBSECTION_FiniteElements,
        "Quantum Graph Spectra" => SUBSECTION_Spectra,
        "spectral_galerkin.md",
        "plots.md"
        ],
    #modules=[MeGraPDE],
    #format = LaTeX(platform = "none"),
    format = Documenter.HTML(prettyurls = false),
    remotes = nothing
)


#= makedocs(;
    root    = "<current-directory>",
    source  = "src",
    build   = "build",
    clean   = true,
    doctest = true,
    modules = Module[],
    repo    = "",
    remotes = "nothing",
    highlightsig = true,
    sitename = "",
    expandfirst = [],
    draft = false
) =#
