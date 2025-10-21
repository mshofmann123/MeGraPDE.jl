# © Anna Weller, University of Cologne, 2023

using LinearAlgebra, Graphs, SparseArrays
using DocStringExtensions

"""
Structure to pass constants used in multigrid iteration across the Level_Parameters.

$(FIELDS)
"""
struct MG_Constants
    # origin vertex of each edge
    origs_e::Vector{Int}
    # terminal vertex of each edge
    terms_e::Vector{Int}
    # Laplacian Matrix
    L::SparseMatrixCSC
    # Degree Matrix
    Deg::SparseMatrixCSC
    # Inverse Degree Matrix
    Deg_inv::SparseMatrixCSC
    # number of edges 
    m::Int
    # number of vertices
    n::Int
end

"""
Structure to pass multigrid parameters.

$(FIELDS)
"""
struct MG_Settings
    # number of pre-smoothing iterations
    nu1::Int
    # number of post-smoothing iterations
    nu2::Int
    # cycle
    mu::Int
end

"""
Structure to pass level discretiation parameters.

$(FIELDS)
"""
struct Level_Parameters
    # metric graph with edge lenghts  ℓ_vec
    ℓ_vec::Vector{Number}
    # vector with inner grid points per edge
    Nx_vec::Vector{Number}
    # vector with step size per edge
    h_vec::Vector{Number}
end


# Intergrid Operators

## nicht-equilateral

"""
    prolongate!(v_fine::Vector, lev_para::Level_Parameters, v_coarse::Vector, mg_const::MG_Constants)

Perform matrix free prolongation of vector 'v_coarse' on level 'Nx_vec/2' to next finer level 'Nx_vec'. 

Calls subroutines 

    prolongate_E!(v_fine::Vector, lev_para::Level_Parameters, v_coarse::Vector, mg_const::MG_Constants)
    prolongate_vector_Nxe(v::Vector ,Nx_e::Int) 

to perform prolongation in block form. 
"""
function prolongate!(v_fine::Vector, lev_para::Level_Parameters, v_coarse::Vector, mg_const::MG_Constants)
    v_fine[1:mg_const.n] = v_coarse[1:mg_const.n] # P_VV * vec_V
    prolongate_E!(v_fine, lev_para, v_coarse, mg_const)   

end
    
"""
    prolongate_E!(v_fine::Vector, lev_para::Level_Parameters, v_coarse::Vector, mg_const::MG_Constants)

Perform prolongation from vertices to edges and inside edges

"""
function prolongate_E!(v_fine::Vector, lev_para::Level_Parameters, v_coarse::Vector, mg_const::MG_Constants)
    count_J = mg_const.n; count_Jm1 = mg_const.n
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j] > 1
            Nxe = Int(lev_para.Nx_vec[j])
            if lev_para.Nx_vec[j]/2 > 1
                v_fine[count_J+1:count_J+Nxe-1] = prolongate_vector_Nxe(v_coarse[count_Jm1+1:count_Jm1+Int(lev_para.Nx_vec[j]/2)-1], Nxe)
            end
            v_fine[count_J+1] += 0.5*v_coarse[mg_const.origs_e[j]]
            v_fine[count_J+Nxe-1] += 0.5*v_coarse[mg_const.terms_e[j]]
            count_J += Nxe-1; count_Jm1 += Int(lev_para.Nx_vec[j]/2)-1;
        end
    end
end

"""
    prolongate_vector_Nxe(v::Vector, Nxe::Int)

Perform prolongation inside edges.

"""
function prolongate_vector_Nxe(v::Vector, Nxe::Int)
    v_out = zeros(Nxe-1);
    v_out[1] = 0.5*v[1];
    count=1;
    for i = 2:2:Nxe-3
        v_out[i] = v[count]
        v_out[i+1] = 0.5*v[count]+0.5*v[count+1]
        count += 1
    end
    v_out[end-1] = v[end]
    v_out[end] = 0.5*v[end]
    return v_out
end
    
"""
    restrict!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)

Perform matrix free restriction of vector 'd_fine' on level 'Nx_vec' to next coarser level 'Nx_vec/2'. 

Calls subroutines 

    restrict_V!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)
    restrict_E!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)
    restrict_vector_Nxe(v::Vector, Nxe::Int)

to perform restriction in block form.

"""
function restrict!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)
    d_coarse[1:mg_const.n] = d_fine[1:mg_const.n] # R_VV * vec_V
    restrict_V!(d_coarse, lev_para, d_fine, mg_const) 
    restrict_E!(d_coarse, lev_para, d_fine, mg_const)
end
    
"""
    restrict_V!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)

Perform restriction from edge to vertex values.

"""
function restrict_V!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)
    counter = mg_const.n+1
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j]<=1
            nothing
        else
            d_coarse[mg_const.origs_e[j]] += 1/2*d_fine[counter]
            d_coarse[mg_const.terms_e[j]] += 1/2*d_fine[counter+Int(lev_para.Nx_vec[j])-2]
            counter += Int(lev_para.Nx_vec[j])-1
        end
    end
end
    
"""
    restrict_E!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)

Perform restrictions inside edges

"""
function restrict_E!(d_coarse::Vector, lev_para::Level_Parameters, d_fine::Vector, mg_const::MG_Constants)
    count_fine=mg_const.n; count_coarse=mg_const.n
    for j = 1:mg_const.m 
        if lev_para.Nx_vec[j]<=1
            nothing
        else
            d_coarse[count_coarse+1:count_coarse+Int(lev_para.Nx_vec[j]/2)-1] = restrict_vector_Nxe(d_fine[count_fine+1:count_fine+Int(lev_para.Nx_vec[j])-1], Int(lev_para.Nx_vec[j]))
            count_fine += Int(lev_para.Nx_vec[j]-1); count_coarse+=Int(lev_para.Nx_vec[j]/2)-1
        end
    end
end
   
"""
    restrict_vector_Nxe(v::Vector,Nxe::Int)

Perform restrictions inside edge

"""
function restrict_vector_Nxe(v::Vector, Nxe::Int)
    v_out = zeros(Int(Nxe/2)-1);
    count = 2
    for i = 1:Int(Nxe/2)-1
        v_out[i] = 0.5*v[count-1]+v[count]+0.5*v[count+1]
        count += 2
    end
    return v_out
end
    
## Integrid (equilateral)
"""
    restrict_vector_e(v::Vector, J::Int)

Equilateral version.
"""
function restrict_vector_e(v::Vector, J::Int)
    v_out = zeros(2^(J-1)-1);
    count = 2
    for i = 1:2^(J-1)-1
        v_out[i] = 0.5*v[count-1]+v[count]+0.5*v[count+1]
        count += 2
    end
    return v_out
end

"""
    prolongate_vector_e(v::Vector, J::Int)

Equilateral version.
"""
function prolongate_vector_e(v::Vector, J::Int)
    v_out = zeros(2^(J)-1);
    v_out[1] = 0.5*v[1]
    count = 1
    for i = 2:2:2^(J)-3
        v_out[i] = v[count]
        v_out[i+1] = 0.5*v[count]+0.5*v[count+1]
        count += 1
    end
    v_out[end-1] = v[end]
    v_out[end] = 0.5*v[end]
    return v_out
end

"""
    restrict_lowest_level!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)

Equilateral version.
"""
function restrict_lowest_level!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)
    d_Jm1[1:mg_const.n] = d_J[1:mg_const.n]
    counter = mg_const.n+1
    for j = 1:mg_const.m
        d_Jm1[mg_const.origs_e[j]] += 1/2*d_J[counter]
        d_Jm1[mg_const.terms_e[j]] += 1/2*d_J[counter+2^J-2]
        counter += 2^J-1
    end
end

"""
    restrict_V!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)

Equilateral version.
"""
function restrict_V!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)
    counter = mg_const.n+1;
    for j = 1:mg_const.m
        d_Jm1[mg_const.origs_e[j]] += 1/2*d_J[counter]
        d_Jm1[mg_const.terms_e[j]] += 1/2*d_J[counter+2^J-2]
        counter += 2^J-1
    end
end

"""
    restrict_E!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)

Equilateral version.
"""
function restrict_E!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)
    count_J = mg_const.n; count_Jm1 = mg_const.n
    for j = 1:mg_const.m 
        d_Jm1[count_J+1:count_J+2^(J-1)-1] = restrict_vector_e(d_J[count_Jm1+1:count_Jm1+2^(J)-1], J)
        count_J += 2^(J-1)-1; count_Jm1 += 2^(J)-1
    end
end

"""
    restrict!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)
    
Equilateral version.
"""
function restrict!(d_Jm1::Vector, J::Int, d_J::Vector, mg_const::MG_Constants)
    # d_Jm1_V
    d_Jm1[1:mg_const.n] = d_J[1:mg_const.n]
    restrict_V!(d_Jm1, J, d_J, mg_const)
    # d_Jm1_E
    restrict_E!(d_Jm1, J, d_J, mg_const)
end

"""
    prolongate_lowest_level!(v_J::Vector, J::Int, v_Jm1::Vector, mg_const::MG_Constants)
    
Equilateral version.
"""
function prolongate_lowest_level!(v_J::Vector, J::Int, v_Jm1::Vector, mg_const::MG_Constants)
    v_J[1:mg_const.n] = v_Jm1[1:mg_const.n]
    count_J = mg_const.n+1; 
    for j = 1:mg_const.m
        v_J[count_J] += 0.5*v_Jm1[mg_const.origs_e[j]]
        v_J[count_J+2^J-2] += 0.5*v_Jm1[mg_const.terms_e[j]]
        count_J+=2^(J)-1; 
    end
end

"""
    prolongate_E!(v_J::Vector, J::Int, v_Jm1::Vector, mg_const::MG_Constants)

Equilateral version.
"""
function prolongate_E!(v_J::Vector, J::Int, v_Jm1::Vector, mg_const::MG_Constants)
    count_J = mg_const.n; count_Jm1 = mg_const.n
    for j=1:mg_const.m
        v_J[count_J+1:count_J+2^(J)-1] = prolongate_vector_e(v_Jm1[count_Jm1+1:count_Jm1+2^(J-1)-1], J)
        v_J[count_J+1] += 0.5*v_Jm1[mg_const.origs_e[j]]
        v_J[count_J+2^J-1] += 0.5*v_Jm1[mg_const.terms_e[j]]
        count_J += 2^(J)-1; count_Jm1 += 2^(J-1)-1; 
    end
end

"""
    prolongate!(v_J::Vector, J::Int, v_Jm1::Vector, mg_const::MG_Constants)

Equilateral version.
"""
function prolongate!(v_J::Vector, J::Int, v_Jm1::Vector, mg_const::MG_Constants)
    v_J[1:mg_const.n] = v_Jm1[1:mg_const.n]
    prolongate_E!(v_J, J, v_Jm1, mg_const)   
end


# Multigrid

## Jacobi iteration

"""
    fill_jacobi_vertices!(u_V::Vector, J::Int, u::Vector, mg_const::MG_Constants)

"""
function fill_jacobi_vertices!(u_V::Vector, J::Int, u::Vector, mg_const::MG_Constants)
    counter = mg_const.n+1
    # only values of first resp. last inner grid points requiered
    for j = 1:mg_const.m
        u_V[mg_const.origs_e[j]] += u[counter]
        u_V[mg_const.terms_e[j]] += u[counter+2^J-2]
        counter += 2^J-1
    end
end

"""
    fill_jacobi_edges!(u::Vector, J::Int, mg_const::MG_Constants)

"""
function fill_jacobi_edges!(u::Vector, J::Int, mg_const::MG_Constants)
    counter = mg_const.n+1;
    for j = 1:mg_const.m
        h = 2.0^(-J);
        scale_Dinv_LowUp = (1/h+2/6*h)^(-1)*(1/6*h-1/h);
        u[counter:counter+2^J-2] = SymTridiagonal(ones(2^J-1),scale_Dinv_LowUp*-0.5*ones((2^J-2)))*u[counter:counter+2^J-2] # hier noch Verbesserungspotenzial --> matrix muss nicht aufgesetellt werden, anwendung von tridiagonalmatrix?!
        u[counter] -= scale_Dinv_LowUp*0.5*u[mg_const.origs_e[j]]
        u[counter+2^J-2] -= scale_Dinv_LowUp*0.5*u[mg_const.terms_e[j]]
        counter += 2^J-1
    end
end

function half_u!(u)
    u = 1/2*u
end

"""
    matrix_free_jacobi!(u::Vector, J::Int, f::Vector, mg_const::MG_Constants)

Iteration of matrix-free Jacobi smoother.     
"""
function matrix_free_jacobi!(u::Vector, J::Int, f::Vector, mg_const::MG_Constants)
    h = 2.0^(-J)
    scale_Dinv = (0.5*(1/h+2/6*h)^(-1))
    scale_Dinv_LowUp = (1/h+2/6*h)^(-1)*(1/6*h-1/h)
    ## apply jacobi to vertex part of u: need to assemble new vector since values of u_V will be needed later on
    u_V = zeros(mg_const.n);                                                                                                           ## work with views of array instead?
    fill_jacobi_vertices!(u_V, J, u, mg_const)
    ## apply jacobi to edge part of u: overwrite u
    fill_jacobi_edges!(u, J, mg_const)
    
    # update vertices values
    u[1:mg_const.n] = u[1:mg_const.n]-scale_Dinv_LowUp*mg_const.Deg_inv*u_V

    # assemble iterate: 1/2 ... u + 1/2 ... f
    u[1:mg_const.n] = half_u!(u[1:mg_const.n])+scale_Dinv*mg_const.Deg_inv*f[1:mg_const.n]
    u[mg_const.n+1:end] = half_u!(u[mg_const.n+1:end])+half_u!(scale_Dinv*(f[mg_const.n+1:end]))
end

"""
    smooth_jacobi!(u0::Vector, J::Int, f_J::Vector, nu1::Int, mg_const::MG_Constants)

Matrix-free Jacobi smoother.     
"""
function smooth_jacobi!(u0::Vector, J::Int, f_J::Vector, nu1::Int, mg_const::MG_Constants)
    for _ = 1:nu1
        matrix_free_jacobi!(u0, J, f_J, mg_const)
    end
end

## Matrix Multiplciations H*u and M*u

"""
    mat_mul_H_VE!(out::Vector, J::Int, u::Vector, mg_const::MG_Constants)

"""
function mat_mul_H_VE!(out::Vector, J::Int, u::Vector, mg_const::MG_Constants)
    counter = mg_const.n+1
    h = 2.0^(-J)
    scale_offdiag = (-1/h+(1/6)*h)
    for j = 1:mg_const.m
        out[mg_const.origs_e[j]] += scale_offdiag*u[counter]
        out[mg_const.terms_e[j]] += scale_offdiag*u[counter+2^J-2]
        counter += 2^J-1
    end
end

"""
    mat_mul_H_EE_EV!(out::Vector, J::Int, u::Vector, mg_const::MG_Constants)

"""
function mat_mul_H_EE_EV!(out::Vector, J::Int, u::Vector, mg_const::MG_Constants)
    counter = mg_const.n+1
    h = 2.0^(-J)
    scale_offdiag = (-1/h+(1/6)*h)
    for j = 1:mg_const.m
        mat_mul_H_EE!(out, J, u, counter)
        # H_EV * u_V
        out[counter] += scale_offdiag*u[mg_const.origs_e[j]]
        out[counter+2^J-2] += scale_offdiag*u[mg_const.terms_e[j]]
        counter += 2^J-1
    end
end

"""
    mat_mul_H_EE!(out::Vector, J::Int, u::Vector, counter::Int)

"""
function mat_mul_H_EE!(out::Vector, J::Int, u::Vector, counter::Int)
    h = 2.0^(-J)
    scale_factor_diag = (1/h+(2/6)*h)
    out[counter:counter+2^(J)-2] = SymTridiagonal(scale_factor_diag*2*ones(2^J-1),(-1/h+1/6*h)*ones(2^J-2))*u[counter:counter+2^(J)-2]
end

"""
    mat_mul_H!(out::Vector, J::Int, u::Vector, mg_const::MG_Constants)

Perform multiplication H_'J'*'u' and store output in 'out'.
"""
function mat_mul_H!(out::Vector, J::Int, u::Vector, mg_const::MG_Constants)
    h = 2.0^(-J)
    scale_diag = (1/h+(2/6)*h)
    # H_VV * u_V
    out[1:mg_const.n] = scale_diag*mg_const.Deg*u[1:mg_const.n]
    # H_VE * u_E
    mat_mul_H_VE!(out, J, u, mg_const)
    # H_EE * u_E
    mat_mul_H_EE_EV!(out, J, u, mg_const)
end


"""
    mat_mul_M!(out::Vector, u::Vector, J::Int, mg_const::MG_Constants)

Perform multiplication M_'J'*'u' and store output in 'out'.
"""
function mat_mul_M!(out::Vector, u::Vector, J::Int, mg_const::MG_Constants)
    h = 2.0^(-J)
    scale_diag = ((2/6)*h)
    # H_VV * u_V
    out[1:mg_const.n] = scale_diag*mg_const.Deg*u[1:mg_const.n]
    # H_VE * u_E
    mat_mul_M_VE!(out, u, J, mg_const)
    # H_EE * u_E
    mat_mul_M_EE_EV!(out, u, J, mg_const)
end

"""
    mat_mul_M_VE!(out::Vector, u::Vector, J::Int, mg_const::MG_Constants)

"""
function mat_mul_M_VE!(out::Vector, u::Vector, J::Int, mg_const::MG_Constants)
    counter = mg_const.n + 1;
    h = 2.0^(-J);
    scale_offdiag = ((1/6)*h);
    for j = 1:mg_const.m
        out[mg_const.origs_e[j]] += scale_offdiag*u[counter]
        out[mg_const.terms_e[j]] += scale_offdiag*u[counter+2^J-2]
        counter += 2^J-1
    end
end

"""
    mat_mul_M_EE_EV!(out::Vector, u::Vector ,J::Int, mg_const::MG_Constants)

"""
function mat_mul_M_EE_EV!(out::Vector, u::Vector ,J::Int, mg_const::MG_Constants)
    counter = mg_const.n+1;
    h = 2.0^(-J);
    scale_offdiag = ((1/6)*h);
    for j = 1:mg_const.m
        mat_mul_M_EE!(out, u, J, counter)
        # H_EV * u_V
        out[counter] += scale_offdiag*u[mg_const.origs_e[j]]
        out[counter+2^J-2] += scale_offdiag*u[mg_const.terms_e[j]]
        counter += 2^J-1
    end
end

"""
    mat_mul_M_EE!(out::Vector, u::Vector, J::Int, counter::Int)

"""
function mat_mul_M_EE!(out::Vector, u::Vector, J::Int, counter::Int)
    h = 2.0^(-J)
    scale_factor_diag = ((2/6)*h)
    out[counter:counter+2^(J)-2] = SymTridiagonal(scale_factor_diag*2*ones(2^J-1),(1/6*h)*ones(2^J-2))*u[counter:counter+2^(J)-2]
end

## MGM

"""
    MGM_matrix_free_jacobi!(u0::Vector, J::Int, f_J::Vector, J0::Int, mg_set::MG_Settings, mg_const::MG_Constants)

Matrix free MGM with jacobi smoother. 
"""
function MGM_matrix_free_jacobi!(u0::Vector, J::Int, f_J::Vector, J0::Int, mg_set::MG_Settings, mg_const::MG_Constants)
    #println(u0)
    # 1. Pre-Smoothing
    smooth_jacobi!(u0, J, f_J, mg_set.nu1, mg_const)
    # 2. Coarse-Grid Correction
    ## 2.1 Residuum
    d_J = zeros(mg_const.n+mg_const.m*(2^J-1));
    mat_mul_H!(d_J, J, u0, mg_const)
    d_J = f_J-d_J
    ## 2.2 Restrict Residuum
    d_Jm1 = zeros(mg_const.n+mg_const.m*(2^(J-1)-1));
    #
    if J==1
        restrict_lowest_level!(d_Jm1, J, d_J, mg_const)
    else
       restrict!(d_Jm1, J, d_J, mg_const)
    end

    ## 2.3 Solve Coarse System

    v_Jm1 = zeros(mg_const.n+mg_const.m*(2^(J-1)-1));
    # solve H_Jm1 v_Jm1 = d_Jm1
    if J-1 == J0  # direct solve
        if J0 === 0
            v_Jm1 = (mg_const.L+abs.(1/6*(mg_const.L+mg_const.Deg)))\d_Jm1
        else
            #FE_1,FE_2=fe_discretization_matrices(G,1,2.0^(-(J-1)))
            #v_Jm1 = Matrix(FE_1+FE_2)\Vector(d_Jm1); 
        end
    else       # μ-times application of MGM_Jm1
        for _ = 1:mg_set.mu
           v_Jm1 = MGM_matrix_free_jacobi!(v_Jm1, J-1, d_Jm1, J0, mg_set, mg_const); 
        end
    end

    ## 2.4 Prolongation
    v_J = zeros(mg_const.n+mg_const.m*(2^(J)-1))
    if J==1
        prolongate_lowest_level!(v_J, J, v_Jm1, mg_const)
    else
        prolongate!(v_J, J, v_Jm1, mg_const)
    end
    # Correction
    u0 = u0+v_J
    # 3. Post-Smoothing
    smooth_jacobi!(u0, J, f_J, mg_set.nu2, mg_const)
    return u0
end

## Solve to discretization error
"""
    solve_mgm(TP::EllipticTestProblem, J::Int; nu1=1, nu2=1, mu=1)

Multigrid solver for elliptic test problem 'TP' at level 'J'. Calls MGM_matrix_free_jacobi!.
"""
function solve_mgm(TP::EllipticTestProblem, J::Int; nu1=1, nu2=1, mu=1)
    # Set up
    mg_set = MG_Settings(nu1, nu2, mu)
    ## assemble graph matrices for direct solver on lowest level (J_0=0)
    L = laplacian_matrix(TP.Γ.G)
    Deg = SparseMatrixCSC(Diagonal(L))
    Deg_inv = SparseMatrixCSC(Diagonal(diag(Deg).^(-1)))
    ## multigrid utils 
    terms_e=zeros(Int,ne(TP.Γ.G)); origs_e=zeros(Int,ne(TP.Γ.G));
    for (m,e) in enumerate(edges(TP.Γ.G))
        terms_e[m] = dst(e);
        origs_e[m] = src(e);
    end
    mg_const = MG_Constants(origs_e, terms_e, L, Deg, Deg_inv, ne(TP.Γ.G), nv(TP.Γ.G))
    #
    h = 2.0^(-J)
    # assemble rhs 
    fe_rhs = zeros(mg_const.n+mg_const.m*(2^J-1))
    f_J = discretize_function(TP.Γ, TP.rhs, Γ.ℓ)
    mat_mul_M!(fe_rhs, f_J, J, mg_const);
    # initialize solution
    u0 = zeros(mg_const.n+mg_const.m*(2^J-1))
    err = 10
    while err > h^2
        u0 = MGM_matrix_free_jacobi!(u0, J, fe_rhs, 0, mg_set, mg_const)
        res = zeros(mg_const.n+mg_const.m*(2^J-1))
        mat_mul_H!(res, J, u0, mg_const)
        res = fe_rhs-res
        err = norm(res,2)
    end
    return u0
end

## Solve with nicmg
## currently works only for solving (L+M)u=f, implement for (L+pot*M)!
"""
    ni_mgm(TP:EllipticTestProblem, J_max::Int)

Nested iteration MG solver for elliptic testproblem 'TP' and Level 'J_max'. Calls MGM_matrix_free_jacobi!.
"""
function ni_mgm(TP::EllipticTestProblem, J_max::Int; nu1=1, nu2=1, mu=1)
    ## currently only implemented for equilateral metric graphs, some changes nec., see cn_mgm
    @assert typeof(TP.Γ) == EquilateralMetricGraph
    ## MG Set up
    mg_set = MG_Settings(nu1, nu2, mu)
    ## assemble graph matrices for direct solver on lowest level (J_0=0)
    L = laplacian_matrix(TP.Γ.G)
    Deg = SparseMatrixCSC(Diagonal(L))
    Deg_inv = SparseMatrixCSC(Diagonal(diag(Deg).^(-1)))
    ## multigrid utils 
    terms_e=zeros(Int,ne(TP.Γ.G)); origs_e=zeros(Int,ne(TP.Γ.G));
    for (m,e) in enumerate(edges(TP.Γ.G))
        terms_e[m] = dst(e);
        origs_e[m] = src(e);
    end
    mg_const = MG_Constants(origs_e, terms_e, L, Deg, Deg_inv, ne(TP.Γ.G), nv(TP.Γ.G))
    ## lowest level
    f_J = discretize_function(TP.Γ, TP.rhs, TP.Γ.ℓ)
    fe_rhs_J = abs.(1/6*(L+Deg))*f_J
    u0 = (L+abs.(1/6*(L+Deg)))\fe_rhs_J
    # prolongate to next level
    u0_Jp1 = zeros(mg_const.n+mg_const.m*(2^1-1))
    prolongate_lowest_level!(u0_Jp1, 0, u0, mg_const)
    # nested iteration
    for J = 1:J_max
        h = 2.0^(-J);
        if J <= 3
            fe_rhs_J = fe_rhs(TP.Γ, TP.rhs, h)
        else
            fe_rhs_J = zeros(mg_const.n+mg_const.m*(2^J-1));
            f_J = discretize_function(TP.Γ, TP.rhs, h)
            mat_mul_M!(fe_rhs_J, f_J, J, mg_const);
        end
        u0 = copy(u0_Jp1);
        err = 10;
        while err > h^2
            u0 = MGM_matrix_free_jacobi!(u0, J, fe_rhs_J, 0, mg_set, mg_const)
            res = zeros(mg_const.n + mg_const.m*(2^J-1))
            mat_mul_H!(res, J, u0, mg_const)
            res = fe_rhs_J-res
            err = norm(res,2)
        end
        # prolongate to next level
        u0_Jp1 = zeros(mg_const.n+mg_const.m*(2^(J+1)-1)) 
        prolongate!(u0_Jp1, J+1, u0, mg_const)
    end
    return u0
end






# CN-Multigrid

## Matrix Multiplication 
"""
    cn_mat_mul!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)

Perform matrix free multiplication (M̂+dt/2*L̂)*'vec'. 

Calls subroutines 

    cn_mat_mul_VV!(out::Vector, vec::Vector, dt::Float64, Nx_vec::Vector, h_vec::Vector, j::Int, counter::Int)
    cn_mat_mul_VE!(out::Vector, vec::Vector,dt::Float64,Nx_vec::Vector,h_vec::Vector,terms_e::Vector,origs_e::Vector,m::Int,counter::Int)
    cn_mat_mul_EV_EE!(out::Vector, vec::Vector,dt::Float64,Nx_vec::Vector,h_vec::Vector,terms_e::Vector,origs_e::Vector,m::Int,counter::Int)
    cn_mat_mul_e!(out::Vector, vec::Vector,dt::Float64,Nx_vec::Vector,h_vec::Vector,j::Int,counter::Int)

to perform the multiplications in block form

    [ (M̂+dt/2*L̂)_VV   (M̂+dt/2*L̂)_VE ]   [ vec_V ]         [(M̂+dt/2*L̂)_VV*vec_V + (M̂+dt/2*L̂)_VE * vec_E ]    
    [                               ]   [       ]   =     [                                            ]
    [ (M̂+dt/2*L̂)_EV   (M̂+dt/2*L̂)_EE ]   [ vec_E ]         [(M̂+dt/2*L̂)_EV*vec_V + (M̂+dt/2*L̂)_EE * vec_E ]

where

                                [ (M̂+dt/2*L̂)_e1                         ]  [ vec_e1 ]
    (M̂+dt/2*L̂)_EE * vec_E =     [                ⋱ ⋱                   ]  [   ⋮    ]  
                                [                        (M̂+dt/2*L̂)_em  ]  [ vec_em ] 

is again performed in block form according to the edges e_1, … ,e_m. 
"""
function cn_mat_mul!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)
    cn_mat_mul_VV!(out, lev_para, vec, dt, mg_const)        # Mat_VV * vec_V
    cn_mat_mul_VE!(out, lev_para, vec, dt, mg_const)        # Mat_VE * vec_E
    cn_mat_mul_EV_EE!(out, lev_para, vec, dt, mg_const)     # Mat_EV * vec_V + Mat_EE * vec_E
end

"""
    cn_mat_mul_VV!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)

Perform multiplication of Mat_VV with vec_V.
"""
function cn_mat_mul_VV!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j] > 1
            h = lev_para.h_vec[j]
            scale_diag = (1/h*(dt/2)+(2/6)*h)
            out[mg_const.origs_e[j]] += scale_diag*vec[mg_const.origs_e[j]]
            out[mg_const.terms_e[j]] += scale_diag*vec[mg_const.terms_e[j]]
        else
            h = lev_para.ℓ_vec[j]
            scale_diag = (1/h*(dt/2)+(2/6)*h)
            out[mg_const.origs_e[j]] += scale_diag*vec[mg_const.origs_e[j]]
            out[mg_const.terms_e[j]] += scale_diag*vec[mg_const.terms_e[j]]
        end
    end
end

"""
    cn_mat_mul_VE!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)

Perform multiplication of Mat_VE with vec_E.
"""
function cn_mat_mul_VE!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)
    counter = mg_const.n+1
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j]>1 
            Nxe = Int(lev_para.Nx_vec[j])
            scale_offdiag = (-1/lev_para.h_vec[j]*(dt/2)+(1/6)*lev_para.h_vec[j])
            out[mg_const.origs_e[j]] += scale_offdiag*vec[counter]
            out[mg_const.terms_e[j]] += scale_offdiag*vec[counter+Nxe-2]
            counter += Nxe-1
        end
    end
end

"""
    cn_mat_mul_EV_EE!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)

Perform multiplication of Mat_EV with vec_V and Mat_EE with vec_E.
"""
function cn_mat_mul_EV_EE!(out::Vector, lev_para::Level_Parameters, vec::Vector, dt::Number, mg_const::MG_Constants)
    counter = mg_const.n+1
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j] > 1
            Nxe = Int(lev_para.Nx_vec[j]); h = lev_para.h_vec[j]
            scale_offdiag = (-1/h*(dt/2)+(1/6)*h)
            cn_mat_mul_e!(out,vec,dt,Nxe,h,counter) # Mat_EE * vec_E is again performed edge wise
            # Mat_EV * u_V
            out[counter] += scale_offdiag*vec[mg_const.origs_e[j]]
            out[counter+Nxe-2] += scale_offdiag*vec[mg_const.terms_e[j]]
            counter += Nxe-1
        end
    end
end

"""
    cn_mat_mul_e!(out::Vector, vec::Vector, dt::Number, Nxe::Int, h::Number, counter::Int)

Perform multiplication of Mat_e with vec_e (for one edge e).
"""
function cn_mat_mul_e!(out::Vector, vec::Vector, dt::Number, Nxe::Int, h::Number, counter::Int)
    scale_factor_diag = (1/h*(dt/2)+(2/6)*h)
    out[counter:counter+Nxe-2] = SymTridiagonal(scale_factor_diag*2*ones(Nxe-1),(-1/h*(dt/2)+1/6*h)*ones(Nxe-2))*vec[counter:counter+Nxe-2]
end

## Jacobi Iteration

"""
    cn_matrix_free_jacobi!(u::Vector, lev_para::Level_Parameters, dt::Number, f::Vector, mg_const::MG_Constants)

Perform weighted jacobi smoother in MGM-CN on 'u' with time step size 'dt', right-hand side 'f' and discretization parameters.
"""
function cn_matrix_free_jacobi!(u::Vector, lev_para::Level_Parameters, dt::Number, f::Vector, mg_const::MG_Constants)

    scale_Dinv = zeros(mg_const.n); 
    cn_weighted_degree!(scale_Dinv, lev_para, dt, mg_const) # assemble weighted degree 

    ## apply jacobi to vertex part of u: need to assemble new vector since values of u_V will be needed later on
    u_V = zeros(mg_const.n);                                                                                                          
    cn_fill_jacobi_vertices!(u_V, lev_para, u, dt, mg_const)
    ## apply jacobi to edge part of u: overwrite u
    cn_fill_jacobi_edges!(u, lev_para, dt, mg_const)
    # update vertex values
    u[1:mg_const.n] = u[1:mg_const.n] - Diagonal(scale_Dinv.^(-1))*u_V
    # assemble iterate: 1/2 ... u + 1/2 ... f
    u[1:mg_const.n] = half_u!(u[1:mg_const.n])+0.5*Diagonal(scale_Dinv.^(-1))*f[1:mg_const.n]
    
    count = mg_const.n+1
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j] > 1
            Nxe = Int(lev_para.Nx_vec[j])
            h = lev_para.ℓ_vec[j]/Nxe
            u[count:count+Nxe-2] = half_u!(u[count:count+Nxe-2])+half_u!(0.5*(((dt/2)*(1/h)+2/6*h)^(-1))*(f[count:count+Nxe-2]))
            count += Nxe-1
        end
    end
end

"""
    cn_weighted_degree!(scale_Dinv::Vector, lev_para:Level_Parameters, dt::Number, mg_const::MG_Constants)

Assemble vector 'scale_Dinv' requiered in cn_matrix_free_jacobi.
"""
function cn_weighted_degree!(scale_Dinv::Vector, lev_para::Level_Parameters, dt::Number, mg_const::MG_Constants)
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j] > 1
            h = lev_para.ℓ_vec[j]/lev_para.Nx_vec[j]
            scale_Dinv[mg_const.origs_e[j]] += (((dt/2)*1/h+2/6*h))
            scale_Dinv[mg_const.terms_e[j]] += (((dt/2)*1/h+2/6*h))
        else
            h = lev_para.ℓ_vec[j]
            scale_Dinv[mg_const.origs_e[j]] += (((dt/2)*1/h+2/6*h))
            scale_Dinv[mg_const.terms_e[j]] += (((dt/2)*1/h+2/6*h))
        end
    end 
end

"""
    cn_fill_jacobi_vertices!(u_V::Vector, lev_para::Level_Parameters, u::Vector, dt::Number, mg_const::MG_Constants)

Perform jacobi smoothing on values on vertices 'u_V'.
"""
function cn_fill_jacobi_vertices!(u_V::Vector, lev_para::Level_Parameters, u::Vector, dt::Number, mg_const::MG_Constants)
    # only values of first resp. last inner grid points requiered
    counter = mg_const.n+1
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j] > 1
            h = lev_para.ℓ_vec[j]/lev_para.Nx_vec[j]
            u_V[mg_const.origs_e[j]] += (1/6*h-(1/h)*(dt/2))*u[counter]
            u_V[mg_const.terms_e[j]] += (1/6*h-(1/h)*(dt/2))*u[counter+Int(lev_para.Nx_vec[j])-2]
            counter += Int(lev_para.Nx_vec[j])-1
        end
    end
end

"""
    cn_fill_jacobi_edges!(u::Vector, lev_para::Level_Parameters, dt::Number, mg_const::MG_Constants)

Perform jacobi smoothing on values on edges.
"""
function cn_fill_jacobi_edges!(u::Vector, lev_para::Level_Parameters, dt::Number, mg_const::MG_Constants)
    counter = mg_const.n+1
    for j = 1:mg_const.m
        if lev_para.Nx_vec[j] > 1
            Nxe = Int(lev_para.Nx_vec[j])
            h = lev_para.ℓ_vec[j]/Nxe;
            scale_Dinv_LowUp = (((1/h)*(dt/2))+2/6*h)^(-1)*(1/6*h-(1/h)*(dt/2));
            u[counter:counter+Nxe-2] = SymTridiagonal(ones(Int(Nxe-1)),scale_Dinv_LowUp*-0.5*ones(Int((Nxe-2))))*u[counter:counter+Nxe-2] # hier noch Verbesserungspotenzial --> matrix muss nicht aufgesetellt werden, anwendung von tridiagonalmatrix?!
            u[counter] -= scale_Dinv_LowUp*0.5*u[mg_const.origs_e[j]]
            u[counter+Nxe-2] -= scale_Dinv_LowUp*0.5*u[mg_const.terms_e[j]]
            counter += Nxe-1
        end
    end
end


## Multigrid

"""
    cn_smooth_jacobi!(u0::Vector, lev_para::Level_Parameters, dt::Number, nu1::Int, f::Vector, mg_const::MG_Constants)

Perform nu1 smooting iterations.
"""
function cn_smooth_jacobi!(u0::Vector, lev_para::Level_Parameters, dt::Number, nu1::Int, f::Vector, mg_const::MG_Constants)
    for _ = 1:nu1
        cn_matrix_free_jacobi!(u0, lev_para, dt, f, mg_const);
    end
end

"""
    cn_coarse_grid_correction_solve!(v_Jm1::Vector, lev_para::Level_Parameters, dt::Float64, d_Jm1::Vector, mg_const::MG_Constants, mg_set::MG_Settings)

"""
function cn_coarse_grid_correction_solve!(v_Jm1::Vector, lev_para::Level_Parameters, dt::Float64, d_Jm1::Vector, mg_const::MG_Constants, mg_set::MG_Settings)
    # solve H_Jm1 v_Jm1 = d_Jm1
    help = lev_para.Nx_vec/2 .<= ones(mg_const.m)
    #
    if help == ones(mg_const.m)  # direct solve
        v_Jm1[1:end] = (((dt/2*mg_const.L+abs.(1/6*(mg_const.L+mg_const.Deg))))\Vector(d_Jm1))[1:end]
    else       # μ-times application of MGM_Jm1 
        for i = 1:mg_set.mu
            lev_para_Jm1 = Level_Parameters(lev_para.ℓ_vec, lev_para.Nx_vec/2, lev_para.h_vec*2)
            cn_mgm_cycle!(v_Jm1, lev_para_Jm1, dt, d_Jm1, mg_const, mg_set)
        end
    end
end

"""
    cn_coarse_grid_correction!(u0::Vector, lev_para::Level_Parameters, dt::Float64, f_J::Vector, mg_const::MG_Constants, mg_set::MG_Settings)

Coarse grid correction including transport to lower level and back.
"""
function cn_coarse_grid_correction!(u0::Vector, lev_para::Level_Parameters, dt::Float64, f_J::Vector, mg_const::MG_Constants, mg_set::MG_Settings)
    # h_vec = lev_para.ℓ_vec./lev_para.Nx_vec;
    ## 2.1 Residuum
    d_J = zeros(length(u0));
    cn_mat_mul!(d_J, lev_para, u0, dt, mg_const);
    d_J = f_J-d_J
    ## 2.2 Restrict Residuum
    d_Jm1 = zeros(mg_const.n+Int(sum([Nxe-1 for Nxe in (lev_para.Nx_vec/2) if Nxe > 1])));
    #
    restrict!(d_Jm1, lev_para, d_J, mg_const)
    ## 2.3 Solve Coarse System
    v_Jm1 = zeros(mg_const.n+Int(sum([Nxe-1 for Nxe in (lev_para.Nx_vec/2) if Nxe > 1])));
    cn_coarse_grid_correction_solve!(v_Jm1, lev_para, dt, d_Jm1, mg_const, mg_set)
    ## 2.4 Prolongation
    v_J = zeros(length(u0))
    prolongate!(v_J, lev_para, v_Jm1, mg_const)
    # Correction
    u0[1:end] = u0[1:end]+v_J[1:end]
end

"""
    cn_mgm_cycle!(u0::Vector, lev_para::Level_Parameters, dt::Float64, f_J::Vector, mg_const::MG_Constants, mg_set::MG_Settings)

Perform one cycle of the CN-MGM method with initial vector 'u0', right-hand side 'f_J', time stepping size 'dt'
"""
function cn_mgm_cycle!(u0::Vector, lev_para::Level_Parameters, dt::Float64, f_J::Vector, mg_const::MG_Constants, mg_set::MG_Settings)
    # 1. Pre-Smoothing
    cn_smooth_jacobi!(u0, lev_para, dt, mg_set.nu1, f_J, mg_const);
    # 2. Coarse-Grid Correction
    cn_coarse_grid_correction!(u0, lev_para, dt, f_J, mg_const, mg_set);
    # 3. Post-Smoothing
    cn_smooth_jacobi!(u0, lev_para, dt, mg_set.nu2, f_J, mg_const);
end

"""
    cn_mgm(TP::TPGeneralizedHeat, T::Number, J::Int; nu1=1, nu2=1, mu=1)

Fully Discretized Scheme: Compute FE-CN-MGM discretization of 'TP' at time 'T' and level 'J'.
"""
function cn_mgm(TP::TPGeneralizedHeat, T::Number, J::Int; nu1=1, nu2=1, mu=1)
    ## MG Set up
    mg_set = MG_Settings(nu1, nu2, mu)
    ## assemble graph matrices for direct solver on lowest level (J_0=0)
    L = laplacian_matrix(TP.Γ.G)
    Deg = SparseMatrixCSC(Diagonal(L))
    Deg_inv = SparseMatrixCSC(Diagonal(diag(Deg).^(-1)))
    ## multigrid utils 
    terms_e=zeros(Int,ne(TP.Γ.G)); origs_e=zeros(Int,ne(TP.Γ.G));
    for (m,e) in enumerate(edges(TP.Γ.G))
        terms_e[m] = dst(e);
        origs_e[m] = src(e);
    end
    mg_const = MG_Constants(origs_e, terms_e, L, Deg, Deg_inv, ne(TP.Γ.G), nv(TP.Γ.G))
    #
    h = 2.0^(-J); dt = h
    # change according to h_max, no uniform level in non-equi case!
    if typeof(TP.Γ) ==EquilateralMetricGraph 
        Nx_vec = TP.Γ.ℓ/h*ones(mg_const.m);
        h_vec = h*ones(mg_const.m);
        ℓ_vec = TP.Γ.ℓ*ones(mg_const.m);
        lev_para = Level_Parameters(ℓ_vec, Nx_vec, h_vec)
        u_t = discretize_function(TP.Γ, TP.u0, h)
    else
        Nx_vec=(2 .^(Int.(ceil.(log.(TP.Γ.ℓ_vec/h)./log.(2)))))
        h_vec = TP.Γ.ℓ_vec./Nx_vec
        lev_para = Level_Parameters(TP.Γ.ℓ_vec, Nx_vec, h_vec)
        u_t = discretize_function(TP.Γ, TP.u0, Nx_vec)
    end
    
    for t = dt:dt:T
        rhs = zeros(length(u_t))
        cn_mat_mul!(rhs, lev_para, u_t, -dt, mg_const)
        res = zeros(length(rhs))
        cn_mat_mul!(res, lev_para, u_t, dt, mg_const)
        while norm(res-rhs) > 0.001*h^2
            cn_mgm_cycle!(u_t, lev_para, dt, rhs, mg_const, mg_set);
            res = zeros(length(rhs))
            cn_mat_mul!(res, lev_para, u_t, dt, mg_const)
        end
    end

    return u_t
end
