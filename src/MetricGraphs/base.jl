# © Anna Weller, University of Cologne, 2023

"""
    metric_graph(G::SimpleGraph, ℓ_vec::Vector; vertex_coords=nothing)

Create metric graph from simple graph 'G' with edge lengths 'ℓ_vec' and optionally assign a coordinate specified in 'coord' to the vertices.

"""
function metric_graph(G::SimpleGraph, ℓ_vec::Vector; vertex_coords=nothing)
    return MetricGraph(G, ℓ_vec, vertex_coords)
end

"""
    metric_graph(G::SimpleGraph, ℓ::Number; vertex_coords=nothing)

Equilateral version with one uniform edge length 'ℓ' assigned to each edge.
"""
function metric_graph(G::SimpleGraph, ℓ::Number; vertex_coords=nothing)
    return EquilateralMetricGraph(G, ℓ, vertex_coords)
end

"""
    edge_length(Γ::MetricGraph, j::Int)

Return edge length of edge 'j.'
"""
function edge_length(Γ::MetricGraph, j::Int)
    return Γ.ℓ_vec[j]
end


"""
    edge_length(Γ::EquilateralMetricGraph)

Equilateral version.
"""
function edge_length(Γ::EquilateralMetricGraph)
    return Γ.ℓ
end


"""
    vol(Γ::MetricGraph)

Return volume ``vol_{\\Gamma} = \\sum_{e \\in \\mathcal{E}} \\ell_e``.
"""
function vol(Γ::MetricGraph)
    return sum(Γ.ℓ_vec)
end

"""
    vol(Γ::EquilateralMetricGraph)

Equilateral version, ``vol = m \\cdot \\ell``.
"""
function vol(Γ::EquilateralMetricGraph)
    return ne(Γ.G)*Γ.ℓ
end
