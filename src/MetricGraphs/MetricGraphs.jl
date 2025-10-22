# © Anna Weller, University of Cologne, 2023

module MetricGraphs

using LinearAlgebra, Graphs, SparseArrays
using DocStringExtensions
using Distributions

import Base: show

export
    # types
    AbstractMetricGraph,
    MetricGraph,
    EquilateralMetricGraph,
    # functions
    metric_graph,
    edge_length,
    vol,
    extended_incidence_matrix,
    extended_metric_graph,
    discretize_function,
    extended_laplacian,
    metric_tree_graph,
    metric_graphene_graph,
    metric_star_graph,
    metric_diamond_graph,
    metric_lollipop_graph,
    metric_barabasi_albert,
    metric_erdos_renyi,
    gaussian_initial,
    model_initial_condition

abstract type AbstractMetricGraph end

function show(io::IO, ::MIME"text/plain", Γ::AbstractMetricGraph)
    if typeof(Γ) == MetricGraph
        return print(io, "{n=$(nv(Γ.G)), m=$(ne(Γ.G))} metric graph")
    elseif typeof(Γ) == EquilateralMetricGraph
        return print(
            io, "{n=$(nv(Γ.G)),m=$(ne(Γ.G)),ℓ=$(edge_length(Γ))} equilateral metric graph"
        )
    end
end

"""
A type representing a Metric Graph (non-equilateral)

$(FIELDS)
"""
struct MetricGraph <: AbstractMetricGraph
    "Simple combinatorial graph"
    G::SimpleGraph
    "Vector containing the edge lengths"
    ℓ_vec::Vector
    "Array containing the coordinates of the vertices; specify 'nothing' if no coordinates available"
    coords::Union{Array,Nothing}

    function MetricGraph(G, ℓ_vec, coords)
        if length(ℓ_vec) != ne(G)
            error("The edge lengths provided does not match the number of edges in G")
        else
            new(G, ℓ_vec, coords)
        end
    end
end

"""
A type representing a Metric Graph with equilateral edge lengths

$(FIELDS)
"""
struct EquilateralMetricGraph <: AbstractMetricGraph
    "Simple combinatorial graph"
    G::SimpleGraph
    "Equilateral edge length"
    ℓ::Number
    "Array containing the coordinates of the vertices; specify 'nothing' if no coordinates available"
    coords::Union{Array,Nothing}
end

function EquilateralMetricGraph(G::SimpleGraph, ℓ::Number)
    return EquilateralMetricGraph(G, ℓ, Nothing)
end

include("base.jl")
include("extended_graph.jl")
include("constructors.jl")

end
