module QuantumGraphSpectra

using LinearAlgebra, Graphs, SparseArrays
using DocStringExtensions
using ..MeGraPDE

import Base:
    show 

export
    # struct
    QuantumGraphDecomposition,
    QuantumGraphEigen,
    EquilateralEigenfunction,
    QuantumGraphEigvals,
    QuantumGraphVertexEigvals,
    EquilateralApproxEigvals,
    EquilateralVertexEigvals,
    EquilateralApproximation,
    # function
    eigen_quantum,
    eigen_quantum_old,
    eigvals_quantum,
    count_eigvals_K,
    eigenfunction,
    equilateral_floor_approximation,
    equilateral_ceil_approximation,
    normalized_laplacian,
    eigvals_equilateral_representation,
    approx_lowest_level,
    nested_iteration_eigenvalue_approximation,
    nested_iteration_newton_trace,
    newton_trace,
    H_matrix


abstract type QuantumGraphDecomposition end 

"""
Spectrum of equilateral graph

$(FIELDS)
"""
struct QuantumGraphEigen <: QuantumGraphDecomposition
    "Number of eigenvalues computed"
    Q::Int
    "eigenvalue"
    Λ::Vector
    "constants A_e"
    A::SparseMatrixCSC
    "constants B_e"
    B::SparseMatrixCSC
end

function show(io::IO, mime::MIME{Symbol("text/plain")}, σ::QuantumGraphDecomposition) 
    #summary(io, σ); println(io)
    println(io, "First $(σ.Q) eigenvalues and eigenfunctions:")
    println(io, "values:")
    show(io, mime, σ.Λ)
    println(io, "\n Coefficients of eigenfunctions ϕ_e = A_e cos(√λ x) + B_e sin (√λ x) for q = 1, …, Q:")
    println(io, "A_e:")
    show(io, mime, σ.A)
    println(io, "\n B_e:")
    show(io, mime, σ.B)
end


"""
Equilateral Eigenfunction

$(FIELDS)
"""
struct QuantumGraphEigfunc
    "eigenvalue"
    λ::Number
    "constants A_e"
    A::Vector
    "constants B_e"
    B::Vector
end

"""
Eigenvalues of quantum graph

$(FIELDS)
"""
struct QuantumGraphEigvals
    "number of eigenvalues computed"
    Q::Int
    "eigenvalues"
    Λ::Vector
end

function show(io::IO, mime::MIME{Symbol("text/plain")}, σ::QuantumGraphEigvals) 
    #summary(io, σ); println(io)
    println(io, "First $(σ.Q) eigenvalues:")
    show(io, mime, σ.Λ)
end

"""
Eigenvalues of quantum graph

$(FIELDS)
"""
struct QuantumGraphVertexEigvals
    "number of eigenvalues computed"
    Q::Int
    "eigenvalues"
    Λ::Vector
end

function show(io::IO, mime::MIME{Symbol("text/plain")}, σ::QuantumGraphVertexEigvals) 
    #summary(io, σ); println(io)
    println(io, "First $(σ.Q) vertex eigenvalues:")
    show(io, mime, σ.Λ)
end


"""
Eigenfunction of quantum graph

$(FIELDS)
"""
struct QuantumGraphEigenfunction
    "eigenvalues"
    λ::Number
    "Coefficients A_e"
    A::Vector
    "Coefficients B_e"
    B::Vector
end

function show(io::IO, mime::MIME{Symbol("text/plain")}, σ::QuantumGraphEigenfunction) 
    #summary(io, σ); println(io)
    println(io, "Eigenfunction corresponding to $(σ.Q) with coefficients:")
    println(io, "A_e:")
    show(io, mime, σ.A)
    println(io, "B_e:")
    show(io, mime, σ.B)
end


"""
EquilateralApproximation of quantum graph

$(FIELDS)
"""
struct EquilateralApproximation
    "Cleaned graph of Equilateral Approximation"
    cl_Γ::MetricGraph
    "Equilateral Approximation"
    Γ̃::EquilateralMetricGraph
    "type of approximation"
    type::Symbol
end



function show(io::IO, mime::MIME{Symbol("text/plain")}, G_h::EquilateralApproximation) 
    #summary(io, σ); println(io)
    println(io, "Equilateral $(G_h.type) approximation:")
    show(io, mime, G_h.Γ̃)
end


"""
Eigenvalues of equilateral approximation

$(FIELDS)
"""
struct EquilateralApproxEigvals
    "Equilateral Approximation"
    Γ̃::EquilateralMetricGraph
    "eigenvalues"
    Λ::Vector
end


include("equilateral_spectrum.jl")
include("non_equilateral_spectrum.jl")

end