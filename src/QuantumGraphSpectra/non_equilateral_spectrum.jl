# © Anna Weller, University of Cologne, 2023

using SparseArrays, LinearAlgebra
using Arpack
using Graphs
using DelimitedFiles

"""
    equilateral_floor_approximation(Γ::MetricGraph, h::Number)

Compute equilateral floor approximation of 'Γ' with equilateral edge length 'h'
"""
function equilateral_floor_approximation(Γ::MetricGraph, h::Number)    
    # -- (cleaned) edge length of approximation --------
    ℓ_approx = zeros(ne(Γ.G))  
    for edge = 1:ne(Γ.G)
        Ne = floor(Γ.ℓ_vec[edge]/h)
        ℓ_approx[edge] = Ne*h
    end
    cl_Γ = MetricGraph(Γ.G,ℓ_approx, nothing)
    # -- incidence matrix of extended graph with equilateraler edge length h ---------
    N_floor = extended_incidence_matrix(MetricGraph(Γ.G, ℓ_approx, nothing), h)
    L_floor = N_floor*N_floor'
    Γ̃ = EquilateralMetricGraph(Graph(L_floor-sparse(Diagonal(L_floor))),h, nothing)
    # ----- Return -----------------------------------------------------------------
    return EquilateralApproximation(cl_Γ, Γ̃, Symbol("floor"))
end

"""
    equilateral_ceil_approximation(Γ::MetricGraph, h::Number)

Compute equilateral ceil approximation of 'Γ' with equilateral edge length 'h'
"""
function equilateral_ceil_approximation(Γ::MetricGraph, h::Number)    
    # -- (cleaned) edge length of approximation --------
    ℓ_approx = zeros(ne(Γ.G))  
    for edge = 1:ne(Γ.G)
        Ne = ceil(Γ.ℓ_vec[edge]/h)
        ℓ_approx[edge] = Ne*h
    end
    cl_Γ = MetricGraph(Γ.G,ℓ_approx, nothing)
    # -- incidence matrix of extended graph with equilateraler edge length h ---------
    N_ceil = extended_incidence_matrix(MetricGraph(Γ.G, ℓ_approx, nothing), h)
    L_ceil = N_ceil*N_ceil'
    Γ̃ = EquilateralMetricGraph(Graph(L_ceil-sparse(Diagonal(L_ceil))),h, nothing)   
    # ----- Return -----------------------------------------------------------------
    return EquilateralApproximation(cl_Γ, Γ̃, Symbol("ceil"))
end

"""
    equilateral_round_approximation(Γ::MetricGraph, h::Number)

Compute rounded equilateral approximation of 'Γ' with equilateral edge length 'h'
"""
function equilateral_round_approximation(Γ::MetricGraph, h::Number)    
    # -- (cleaned) edge length of approximation --------
    ℓ_approx = zeros(ne(Γ.G))  
    for edge = 1:ne(Γ.G)
        Ne = round(Γ.ℓ_vec[edge]/h)
        ℓ_approx[edge] = Ne*h
    end
    cl_Γ = MetricGraph(Γ.G,ℓ_approx, nothing)
    # -- incidence matrix of extended graph with equilateraler edge length h ---------
    N_ceil = extended_incidence_matrix(MetricGraph(Γ.G, ℓ_approx, nothing), h)
    Γ̃ = EquilateralMetricGraph(Graph(N_ceil*N_ceil'),h, nothing)
    # ----- Return -----------------------------------------------------------------
    return EquilateralApproximation(cl_Γ, Γ̃, Symbol("round"))
end



function normalized_laplacian(Inc::SparseMatrixCSC)
    L = Inc*Inc';
    L_norm = sparse(Diagonal(diag(L).^(-1/2)))*L*sparse(Diagonal(diag(L).^(-1/2)));
    return L_norm
end


function transform_eigvals(Γ::EquilateralMetricGraph, μ_vec::Vector, k::Int)
    return ((1/Γ.ℓ)*acos.(ones(length(μ_vec))-μ_vec)+k*pi*ones(length(μ_vec))).^2
end

function transform_eigvals(Γ::EquilateralMetricGraph, μ::Number, k::Int)
    return ((1/Γ.ℓ)*(acos(1-μ)+k*pi))^2
end


"""
    eigvals_equilateral_representation(Γ::MetricGraph, h::Number)

Compute the exact eigenvalues of 'Γ' by an equilateral representation with 
edge length 'h'
"""
function eigvals_equilateral_representation(Γ::MetricGraph, h::Number)
    if isinteger.(Γ.ℓ_vec/h) != BitVector(ones(ne(Γ.G)))
        error("$h is not a comman divisor of all edge lengths")
    else
        Inc = extended_incidence_matrix(Γ, h);
        L_norm = normalized_laplacian(Inc);
        μ_vec = eigvals(Matrix(L_norm));
        #
        Γ̃ = EquilateralMetricGraph(Graph(Inc*Inc'), h, nothing)
        if isapprox(μ_vec[end],2)
            return QuantumGraphEigvals(length(μ_vec[2:end-2]), transform_eigvals(Γ̃, μ_vec[2:end-2], 0))
        else
            return QuantumGraphEigvals(length(μ_vec[2:end-1]), transform_eigvals(Γ̃, μ_vec[2:end-1], 0))
        end
    end
end


"""
    approx_lowest_level(Γ::MetricGraph, h_min::Number; Q=2)

Compute eigenvalue approximations by equilateral ceil and floor approximations of the first 'Q' eigenvalues at the lowest discretization level 'h_min' in the nested iteration.
""" 
function approx_lowest_level(Γ::MetricGraph, h_min::Number, Q::Int)
    Γ_floor = equilateral_floor_approximation(Γ, h_min) 
    L_floor_norm = normalized_laplacian(incidence_matrix(Γ_floor.Γ̃.G, oriented=true))
    μ_floor = eigs(L_floor_norm, nev=Q, maxiter=100000, which=:SR)[1][2:Q];
    Γ_ceil = equilateral_ceil_approximation(Γ, h_min) 
    L_ceil_norm = normalized_laplacian(incidence_matrix(Γ_floor.Γ̃.G, oriented=true))
    μ_ceil = eigs(L_ceil_norm, nev=Q, maxiter=100000, which=:SR)[1][2:Q];
    return transform_eigvals(Γ_floor.Γ̃, μ_floor, 0), transform_eigvals(Γ_ceil.Γ̃, μ_ceil, 0)
end

"""
    nested_iteration_eigenvalue_approximation(Γ::MetricGraph; lev_zero=0, lev_max=7, Q=2, save_each_lev=false)

Approximate first 'Q' eigenvalues of 'Γ' via equilateral approximations using a nested itertation approach.
"""
function nested_iteration_eigenvalue_approximation(Γ::MetricGraph, lev_zero, lev_max, Q, save_each_lev)
    h_min = 2.0^(-lev_zero)
    Λ_floor,Λ_ceil = approx_lowest_level(Γ,h_min,Q)
    #
    if save_each_lev==true
        writedlm("lambda_floor_lev_zero.txt", Λ_floor)
        writedlm("lambda_ceil_lev_zero.txt", Λ_ceil)
    end
    #
    for lev = lev_zero+1:lev_max
        h = 2.0^(-lev)
        Λ_floor_lev=[]; Λ_ceil_lev=[];
        for q = 1:Q-1
            # initial guess for sparse eigenvalue solver
            initial_guess = ((1-cos(sqrt(Λ_floor[q])*h))+(1-cos(sqrt(Λ_ceil[q])*h)))/2
            # sparse eigenvalue solver
            Γ_floor = equilateral_floor_approximation(Γ, h)   
            L_floor_norm = normalized_laplacian(incidence_matrix(Γ_floor.Γ̃.G, oriented=true))
            μ_floor_q = eigs(L_floor_norm,nev=3,sigma=initial_guess+10^(-8))[1][1:3]
            for i = 1:3
                push!(Λ_floor_lev, ((1/h)*acos(1-μ_floor_q[i]))^2)
            end
            Γ_ceil = equilateral_ceil_approximation(Γ, h)   
            L_ceil_norm = normalized_laplacian(incidence_matrix(Γ_ceil.Γ̃.G, oriented=true))
            μ_ceil_q = eigs(L_ceil_norm,nev=3,sigma=initial_guess+10^(-8))[1][1:3]
            for i = 1:3
                push!(Λ_ceil_lev, ((1/h)*acos(1-μ_ceil_q[i]))^2)
            end
        end

        Λ_floor[1:Q-1] = sort((unique(trunc.(Λ_floor_lev, digits = 10))))[1:Q-1]
        Λ_ceil[1:Q-1] = sort((unique(trunc.(Λ_ceil_lev, digits = 10))))[1:Q-1]
        
        if save_each_lev==true
            writedlm("lambda_floor_lev_$lev.txt", Λ_floor)
            writedlm("lambda_ceil_lev_$lev.txt", Λ_ceil)
        end
    end

    return Λ_floor, Λ_ceil
end


### NEP

"""
    H_matrix(z::Number, Γ::MetricGraph)

Compute H(z) for a metric graph 'Γ'.
"""
function H_matrix(z::Number, Γ::MetricGraph)
    bfN = incidence_matrix(Γ.G, oriented=true)
    # Diagonal Entries: Weighted degree matrix with edge weights -cot(√z * ℓ_e)
    f_diag= z -> -cot.(sqrt(z)*Γ.ℓ_vec);
    W1 = Diagonal(f_diag(z));
    H_diag = Diagonal(diag(bfN*W1*bfN'));

    # Off-Diagonal Entries: Weighted Adjacency matrix with edge weights 1/sin(√ z * ℓ_e)
    f_offdiag = z -> 1 ./ sin.(sqrt(z)*Γ.ℓ_vec); 
    W2 = Diagonal(f_offdiag(z));
    H_off_diag = -(bfN*W2*bfN'-Diagonal(diag(bfN*W2*bfN'))) ## inner bracket gives off-diagonal entries of weighted Laplacian, this is - Adjacency

    return H_diag + H_off_diag 
end

"""
    H_matrix(z::Number, bfN::SparseMatrixCSC, ℓ_vec::Vector)

Compute H(z) for a graph with incidence matrix 'Inc' and edge length 'ℓ_vec'.
"""
function H_matrix(z::Number, bfN::SparseMatrixCSC, ℓ_vec::Vector)
    # Diagonal Entries: Weighted degree matrix with edge weights -cot(√z * ℓ_e)
    f_diag= z -> -cot.(sqrt(z)*ℓ_vec);
    W1 = Diagonal(f_diag(z));
    H_diag = Diagonal(diag(bfN*W1*bfN'));

    # Off-Diagonal Entries: Weighted Adjacency matrix with edge weights 1/sin(√ z * ℓ_e)
    f_offdiag = z -> 1 ./ sin.(sqrt(z)*ℓ_vec); 
    W2 = Diagonal(f_offdiag(z));
    H_off_diag = -(bfN*W2*bfN'-Diagonal(diag(bfN*W2*bfN'))) ## inner bracket gives off-diagonal entries of weighted Laplacian, this is - Adjacency

    return H_diag + H_off_diag 
end

"""
    H_matrix_deriv(z::Number, Γ::MetricGraph)

Compute H'(z) for a metric graph with incidence matrix 'Inc' and edge length 'ℓ_vec'.

"""
function H_matrix_deriv(z::Number, Γ::MetricGraph)
    bfN = incidence_matrix(Γ.G, oriented = true)
    df_diag = z -> Γ.ℓ_vec ./(2*sqrt(z)*sin.(sqrt(z).*Γ.ℓ_vec).^2)
    dW1 = Diagonal(df_diag(z))
    H_deriv_diag = Diagonal(diag(bfN*dW1*bfN'))

    df_offdiag = z -> -(Γ.ℓ_vec.*cot.(sqrt(z).*Γ.ℓ_vec)./(2*sqrt(z)*sin.(sqrt(z).*Γ.ℓ_vec)))
    dW2 = Diagonal(df_offdiag(z))
    H_deriv_off_diag = -(bfN*dW2*bfN'-Diagonal(diag(bfN*dW2*bfN')))

    return H_deriv_diag+H_deriv_off_diag
end

"""
    H_matrix_deriv(z::Number, bfN::SparseMatrixCSC, ℓ_vec::Vector)

Compute H'(z) for a metric graph with incidence matrix 'Inc' and edge length 'ℓ_vec'.

"""
function H_matrix_deriv(z::Number, bfN::SparseMatrixCSC, ℓ_vec::Vector)
    df_diag = z -> ℓ_vec ./(2*sqrt(z)*sin.(sqrt(z).*ℓ_vec).^2)
    dW1 = Diagonal(df_diag(z))
    H_deriv_diag = Diagonal(diag(bfN*dW1*bfN'))

    df_offdiag = z -> -(ℓ_vec.*cot.(sqrt(z).*ℓ_vec)./(2*sqrt(z)*sin.(sqrt(z).*ℓ_vec)))
    dW2 = Diagonal(df_offdiag(z))
    H_deriv_off_diag = -(bfN*dW2*bfN'-Diagonal(diag(bfN*dW2*bfN')))

    return H_deriv_diag+H_deriv_off_diag
end


### Newton-trace Iteration 

"""
    newton_trace(Γ::MetricGraph, z_start::Number)

Newton-trace iteration to determine roots of det(H(z)).

"""
function newton_trace(Γ::MetricGraph, z_start::Number)
    z = z_start
    tol = 1e-10; maxit = 1000; 
    iter = 0;
    while iter <= maxit && 1/cond(H_matrix(z,Γ),1) > tol
        corr = tr(Matrix(H_matrix(z,Γ))\H_matrix_deriv(z,Γ))
        z = z-1/corr
        iter += 1
    end
    # end
    return z, iter
end

"""
    norm_equilateral_eigenfunction(ϕ::EquilateralEigfunc)

Compute norm of eigenfunction 'ϕ' of equilateral metric graph.
"""
function eigfunc_norm(λ::Number, A::Vector, B::Vector, ℓ_vec::Vector, m::Int)
    norm = 0
    λ_sqrt = sqrt(λ)
    for j = 1:m
        Ae = A[j]; Be = B[j]
        norm += (2*λ_sqrt*ℓ_vec[j]*(Ae^2+Be^2)+(Ae^2-Be^2)*sin(2*λ_sqrt*ℓ_vec[j])+2*Ae*Be*(1-cos(2*λ_sqrt*ℓ_vec[j])))/(4*λ_sqrt)
    end
    return sqrt(norm)
end

"""
    nested_iteration_newton_trace(Γ::MetricGraph; lev_zero=0, lev_max=7, Q=5, save_each_lev=false, return_eigvecs=false)

Conduct nested iteration newton trace algorithm to find the first 'Q' eigenvalues of 'Γ'.

"""
function nested_iteration_newton_trace(Γ::MetricGraph; lev_zero=0, lev_max=7, Q=5, save_each_lev=false, return_eigvecs=false)
    λ_floor,λ_ceil = nested_iteration_eigenvalue_approximation(Γ, lev_zero, lev_max, Q, save_each_lev)
    Λ=[]; iterations=[]
    for q = 1:Q-1
        z, iter = newton_trace(Γ, (λ_floor[q]+λ_ceil[q])/2)
        push!(Λ,z); push!(iterations,iter)
    end
    if return_eigvecs == false
        return Λ
    else 
        A = spzeros(ne(Γ.G),Q-1); B = spzeros(ne(Γ.G),Q-1)
        for q=1:Q-1
            λ = Λ[q]
            Φ = eigs(H_matrix(Λ[q], Γ),sigma=0,nev=1)[2]
            A_q = zeros(ne(Γ.G)); B_q = zeros(ne(Γ.G))
            for (j,e) in enumerate(edges(Γ.G))
                e_0 = src(e); e_ℓ = dst(e) 
                A_q[j] = Φ[e_0]
                B_q[j] = 1/(sin(sqrt(λ)*Γ.ℓ_vec[j]))*(Φ[e_ℓ]-Φ[e_0]*cos(sqrt(λ)*Γ.ℓ_vec[j]))
            end
            norm = eigfunc_norm(λ,A_q,B_q,Γ.ℓ_vec,ne(Γ.G))
            A[:,q] = A_q./norm; B[:,q] = B_q./norm
        end
        return QuantumGraphEigen(Q-1, Λ, A, B)
    end
end




function eigvals_quantum(Γ::MetricGraph)
    error("eigvals_quantum is currently only implemented for equilateral metric graphs. 
    Please decide according to you situation which routine is best suited:
        - eigvals_equilateral_representation(Γ::MetricGraph, h::Number) can be used if Γ can be represented by an equilateral representation with a moderate number of vertices
        - equilateral_approximation_newton_trace(Γ::MetricGraph; lev=5) can be used for small graphs 
        - nested_iteration_newton_trace(Γ::MetricGraph, lev_zero=0, lev_max=7, Q=2) can be applied for large graphs with large deviations in the edge lengths to find the first 'Q' eigenvalues
    ")
end