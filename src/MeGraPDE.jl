# © Anna Weller, University of Cologne, 2023

module MeGraPDE

using DocStringExtensions
using SparseArrays

import Base:
    show 

export 
    # MetricGraphs
    MetricGraph,
    EquilateralMetricGraph,
    metric_graph,
    edge_length,
    vol,
    extended_incidence_matrix,
    extended_laplacian,
    discretize_function,
    metric_tree_graph,
    metric_graphene_graph,
    metric_star_graph,
    metric_diamond_graph,
    metric_lollipop_graph,
    metric_barabasi_albert,
    metric_erdos_renyi,
    gaussian_initial,
    model_initial_condition,
    # finite elements 
    fe_matrices,
    fe_rhs,
    fe_discretization,
    finite_element_solver,
    fe_error,
    # multigrid 
    cn_mgm_iter!,
    cn_mgm,
    # spectra equilateral
    QuantumGraphEigen,
    eigen_quantum,
    eigen_quantum_old,
    eigvals_quantum,
    count_eigvals_K,
    eigenfunction,
    # spectra non-equilateral
    equilateral_floor_approximation,
    equilateral_ceil_approximation,
    normalized_laplacian,
    eigvals_equilateral_representation,
    approx_lowest_level,
    nested_iteration_eigenvalue_approximation,
    nested_iteration_newton_trace,
    newton_trace,
    H_matrix,
    # spectral galerkin
    projection_coefficients,
    projection_coefficients_filon,
    spectral_solution,
    spectral_solver,
    L2_norm_spectral,
    H1_seminorm_spectral,
    # Testproblems
    EllipticTestProblem,
    ParabolicTestProblem,
    TPGeneralizedHeat,
    TPReactionDiffusion,
    TPPoisson,
    TestProblem242,
    TestProblem243,
    TestProblem244,
    TestProblem245,
    TestProblem711,
    TestProblem721,
    test_problem_6_1_6,
    # Plots
    plot_graph_3d,
    plot_function_3d,
    # animations
    animate_diffusion,
    animate_fractional_diffusion

    
# Metric Graphs
include("MetricGraphs/MetricGraphs.jl")
using .MetricGraphs

# Spectra
include("QuantumGraphSpectra/QuantumGraphSpectra.jl")
using .QuantumGraphSpectra

# Testproblems 
include("testproblems.jl")

# Finite Elements
include("finite_elements/finite_element_discretization.jl")
include("finite_elements/multigrid.jl")
include("finite_elements/multigrid_matrix_free.jl")

# Spectral Galerkin
include("spectral_galerkin/projection_coefficients.jl")
include("spectral_galerkin/spectral_solution.jl")

# plots
include("plots/graph_plots.jl")
include("examples/animations.jl")

end
