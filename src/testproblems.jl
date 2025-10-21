# © Anna Weller, University of Cologne, 2023a

using Graphs

# Testproblems & Examples

abstract type TestProblem end

function show(io::IO, mime::MIME{Symbol("text/plain")}, TP::TestProblem)
    println(io, "Test Problem on Γ")
    println(io, "Γ:")
    show(io, mime, TP.Γ)
end

"""
Elliptic Test Problem 

$(FIELDS)
"""
struct EllipticTestProblem <: TestProblem
    "Potential"
    pot::Number
    "Metric Graph"
    Γ::Union{MetricGraph,EquilateralMetricGraph}
    "Right-hand side"
    rhs::Vector{Function}
    "Exact solution"
    u::Vector{Function}
    "Derivative of exact solution"
    u_deriv::Vector{Function}
end

"""
    test_problem_2_4_2()

Elliptic test problem on 5-star graph with equilateral edge length π+π/2.
"""
function test_problem_2_4_2()
    Γ = metric_star_graph(; ℓ=π+π/2)

    # right hand sinde
    rhs = [x -> -6*sin(x), x -> 2*sin(x), x -> 2*sin(x), x -> 2*sin(x)];

    # exact solution
    u = [x -> -3*sin(x), x -> sin(x), x -> sin(x), x -> sin(x)];

    # derivative of exact solution (for H^1 norm)
    u_deriv = [x -> -3*cos(x), x -> cos(x), x -> cos(x), x -> cos(x)];

    return EllipticTestProblem(1, Γ, rhs, u, u_deriv)
end

"""
    TestProblem242

Elliptic test problem on 5-star graph with equilateral edge length π+π/2.
"""
const TestProblem242 = test_problem_2_4_2()

"""
    test_problem_2_4_3()

Elliptic test problem diamond graph with equilateral edge length 2π.
"""
function test_problem_2_4_3()
    # diamond graph
    Γ = metric_diamond_graph(; ℓ=2*π)
    # right hand sinde
    rhs = [x -> 2*sin(x), x -> -4*sin(x), x -> 2*sin(x), x -> 2*sin(x), x -> -2*sin(x)];

    # exact solution
    u = [x -> sin(x), x -> -2*sin(x), x -> sin(x), x -> sin(x), x -> -sin(x)];

    # derivative of exact solution (for H^1 norm)
    u_deriv = [x -> cos(x), x -> -2*cos(x), x -> cos(x), x -> cos(x), x -> -cos(x)];

    return EllipticTestProblem(1, Γ, rhs, u, u_deriv)
end

"""
    TestProblem243  

Elliptic test problem diamond graph with equilateral edge length 2π.
"""
const TestProblem243 = test_problem_2_4_3()

"""
    test_problem_7_1_1()

```math
\\mathcal{H} u + u = f
```
on star graph with
```math
f = -\\frac{(\\exp(-(x-x_0)^2/s^2) \\cdot (4 (x-x_0)^2-2s^2)}{s^4} + \\exp \\left(- \\frac{(x-x_0)^2}{s^2} \\right)
```
on one edge and exact solution
```math
\\exp \\left(- \\frac{(x-x_0)^2}{s^2} \\right).
```
"""
function test_problem_7_1_1()
    Γ = metric_star_graph()
    # choose edge to assign non-zero function on
    x0 = Γ.ℓ/2;
    s=Γ.ℓ/10;
    edge_init = rand(1:ne(Γ.G))
    # edge_init = 1
    # assemble vectors of functions
    rhs = Vector{Function}(undef, ne(Γ.G));
    u = Vector{Function}(undef, ne(Γ.G));
    u_deriv = Vector{Function}(undef, ne(Γ.G));
    # assign functions
    for j in 1:ne(Γ.G)
        if j == edge_init
            rhs[j] = x -> -(exp(-(x-x0)^2/s^2)*(4*(x-x0)^2-2*s^2))/s^4+exp(-(x-x0)^2/s^2);
            u[j] = x -> exp(-(x-x0)^2/s^2);
            u_deriv[j] = x -> -(2*(x-x0)*exp(-(x-x0)^2/s^2))/s^2
        else
            rhs[j] = x -> 0;
            u[j] = x -> 0;
            u_deriv[j] = x -> 0;
        end
    end
    return EllipticTestProblem(1, Γ, rhs, u, u_deriv)
end

"""
    TestProblem711

Elliptic test problem on star graph.
"""
const TestProblem711 = test_problem_7_1_1()

## ToDo: Write function list_testproblems() that lists all test problems with short description. Same for example graphs?

## ToDo: Write function plot_testproblem() that illustrates test problems --> elliptic print rhs?, parabolic print init

"""
Parabolic Test Problem 

$(FIELDS)
"""
struct ParabolicTestProblem
    "Metric Graph"
    Γ::Union{MetricGraph,EquilateralMetricGraph}
    "Initial Condition"
    u0::Vector{Function}
    "Right-hand side"
    rhs::Union{Vector{Function},Nothing}
    "Reaction term"
    react::Union{Vector{Function},Nothing}
    "Exact solution"
    u::Union{Vector{Function},Nothing}
    "Derivative of exact solution"
    u_deriv::Union{Vector{Function},Nothing}
end

"""
Test Problem related to Generalized Heat Equation on Γ.

$(FIELDS)
"""
struct TPGeneralizedHeat <: TestProblem
    "Metric Graph"
    Γ::Union{MetricGraph,EquilateralMetricGraph}
    "Initial Condition"
    u0::Vector{Function}
    "right-hand side"
    rhs::Union{Vector{Function},Nothing}
    "Exact solution"
    u::Union{Vector{Function},Nothing}
    "Derivative of exact solution"
    u_deriv::Union{Vector{Function},Nothing}
end

"""
    test_problem_2_4_4()

Parabolic Test Problem on equilateral star graph with eigenfunction as initial condition.
"""
function test_problem_2_4_4(; eigfunc_no=5)
    Γ=metric_star_graph()
    # compute spectrum for initial condition
    σ = eigen_quantum(Γ; K=2)
    A_q = σ.A[:, eigfunc_no];
    B_q = σ.B[:, eigfunc_no];
    λ_q = σ.Λ[eigfunc_no]
    # initial condition & exact solution
    u0 = [];
    u = [];
    u_deriv = []
    for j in 1:ne(Γ.G)
        push!(u0, x -> A_q[j]*cos(sqrt(λ_q)*x)+B_q[j]*sin(sqrt(λ_q)*x));
        push!(u, (t, x) -> exp(-t*λ_q)*(A_q[j]*cos(sqrt(λ_q)*x)+B_q[j]*sin(sqrt(λ_q)*x)));
        push!(
            u_deriv,
            (t, x) ->
                exp(
                    -t*λ_q
                )*(-A_q[j]*sqrt(λ_q)*sin(sqrt(λ_q)*x)+B_q[j]*sqrt(λ_q)*cos(sqrt(λ_q)*x)),
        );
    end
    # return
    return TPGeneralizedHeat(Γ, u0, nothing, u, u_deriv)
end

"""
Heat equation on star graph with eigenfunction ``\\phi_5`` as initial condition.
"""
const TestProblem244 = test_problem_2_4_4()

"""
    test_problem_2_4_5()

Parabolic Test Problem on equilateral diamond graph with eigenfunction as initial condition.
"""
function test_problem_2_4_5(; eigfunc_no=3)
    Γ = metric_diamond_graph()
    # compute spectrum for initial condition
    σ = eigen_quantum(Γ; K=2)
    A_q = σ.A[:, eigfunc_no];
    B_q = σ.B[:, eigfunc_no];
    λ_q = σ.Λ[eigfunc_no]
    # initial condition & exact solution
    u0 = [];
    u = [];
    u_deriv = []
    for j in 1:ne(Γ.G)
        push!(u0, x -> A_q[j]*cos(sqrt(λ_q)*x)+B_q[j]*sin(sqrt(λ_q)*x));
        push!(u, (t, x) -> exp(-t*λ_q)*(A_q[j]*cos(sqrt(λ_q)*x)+B_q[j]*sin(sqrt(λ_q)*x)));
        push!(
            u_deriv,
            (t, x) ->
                exp(
                    -t*λ_q
                )*(-A_q[j]*sqrt(λ_q)*sin(sqrt(λ_q)*x)+B_q[j]*sqrt(λ_q)*cos(sqrt(λ_q)*x)),
        );
    end
    # return
    return TPGeneralizedHeat(Γ, u0, nothing, u, u_deriv)
end

"""
Heat equation on diamond graph with eigenfunction ``\\phi_3`` as initial condition.
"""
const TestProblem245 = test_problem_2_4_5()

function test_problem_6_1_6(n::Int, k::Int)
    Γ = metric_barabasi_albert(n, k; ℓ=:non_equi)
    σ = nested_iteration_newton_trace(Γ; lev_zero=0, lev_max=3, Q=2, return_eigvecs=true)
    A_q = σ.A[:, 1];
    B_q = σ.B[:, 1];
    λ_q = σ.Λ[1]
    # initial condition & exact solution
    u0 = [];
    u = [];
    u_deriv = []
    for j in 1:ne(Γ.G)
        push!(u0, x -> A_q[j]*cos(sqrt(λ_q)*x)+B_q[j]*sin(sqrt(λ_q)*x));
        push!(u, (t, x) -> exp(-t*λ_q)*(A_q[j]*cos(sqrt(λ_q)*x)+B_q[j]*sin(sqrt(λ_q)*x)));
        push!(
            u_deriv,
            (t, x) ->
                exp(
                    -t*λ_q
                )*(-A_q[j]*sqrt(λ_q)*sin(sqrt(λ_q)*x)+B_q[j]*sqrt(λ_q)*cos(sqrt(λ_q)*x)),
        );
    end
    # return
    return TPGeneralizedHeat(Γ, u0, nothing, u, u_deriv)
end

"""
    test_problem_7_2_1()    

"""
function test_problem_7_2_1()
    Γ = metric_star_graph()
    x0=Γ.ℓ/2;
    s=Γ.ℓ/10;
    edge_init = rand(1:ne(Γ.G))
    # edge_init = 1
    # assemble vectors of functions
    u0 = Vector{Function}(undef, ne(Γ.G));
    # assign functions
    for j in 1:ne(Γ.G)
        if j == edge_init
            u0[j] = x -> exp(-(x-x0)^2/s^2);
        else
            #u0[j] = x -> exp(-(x-x0)^2/s^2);
            u0[j] = x -> 0;
        end
    end
    return TPGeneralizedHeat(Γ, u0, nothing, nothing, nothing)
end

"""
Heat equation with Gaussian-type inital condition.
"""
const TestProblem721 = test_problem_7_2_1()

"""
Test Problem related to Reaction-Diffusion Equation on Γ.

$(FIELDS)
"""
struct TPReactionDiffusion <: TestProblem
    "Metric Graph"
    Γ::Union{MetricGraph,EquilateralMetricGraph}
    "Initial Condition"
    u0::Vector{Function}
    "reaction term"
    react::Union{Vector{Function},Nothing}
    "Exact solution"
    u::Union{Vector{Function},Nothing}
    "Derivative of exact solution"
    u_deriv::Union{Vector{Function},Nothing}
end
