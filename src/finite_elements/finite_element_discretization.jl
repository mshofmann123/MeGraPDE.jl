# © Anna Weller, University of Cologne, 2023

using  BSplines

function fe_discretization_parameter(Γ::AbstractMetricGraph, h_max::Number)
    if typeof(Γ) == EquilateralMetricGraph
        return Int(ceil(Γ.ℓ/h_max));
    elseif typeof(Γ) == MetricGraph
        return Int.(ceil(Γ.ℓ_vec/h_max));
    end
end

"""
    fe_matrices(Γ::MetricGraph, h_max::Number; mass_approx = false)

Assemble finite element stiffness and mass matrix.

Return (stiffness matrix, mass matrix). Set mass_approx = true to use mass matrix approximation via Trapezoidal rule.
"""
function fe_matrices(Γ::MetricGraph, h_max::Number; mass_approx = false)
    # initialize weight matrices 
    W = spzeros(0,0); W_inv = spzeros(0,0);
    for j=1:ne(Γ.G);
        Nx_e = Int.(ceil(Γ.ℓ_vec[j]/h_max));
        h_e = Γ.ℓ_vec[j]/Nx_e;
        W = blockdiag(W,1/h_e*sparse(I(Nx_e)));
        W_inv = blockdiag(W_inv,h_e*sparse(I(Nx_e)));
    end
    Inc = extended_incidence_matrix(Γ,h_max);
    if mass_approx == false
        # stiffness matrix  L = Inc*W*Inc'
        # mass matrix       M = 1/6*(abs.(Inc*W_inv*Inc')+abs.(Diagonal(Inc*W_inv*Inc')))
        # return L,M
        return SparseMatrixCSC(Inc*W*Inc'),SparseMatrixCSC(1/6*(abs.(Inc*W_inv*Inc')+abs.(Diagonal(Inc*W_inv*Inc'))))
    else
        # stiffness matrix  L = Inc*W*Inc'
        # mass matrix       M = Diagonal(1/2*diag(Inc*W_inv*Inc'))
        # return L,M 
        return SparseMatrixCSC(Inc*W*Inc'),Diagonal(1/2*diag(Inc*W_inv*Inc'))
    end
end


""" 
    fe_matrices(Γ::EquilateralMetricGraph, h_max::Number; mass_approx = false)

When called for equilateral graphs, apply uniform discretization on each edge.

"""
function fe_matrices(Γ::EquilateralMetricGraph, h_max::Number; mass_approx = false)
    # number of grid points on each edge
    Nx = Int(ceil(Γ.ℓ/h_max));
    h = Γ.ℓ/Nx;
    # weight matrix
    W = Diagonal(1/h*ones(ne(Γ.G)*Nx)) 
    W_inv = Diagonal(h*ones(ne(Γ.G)*Nx)) 
    # incidencte matrix
    Inc = extended_incidence_matrix(Γ,h);
    if mass_approx == false
        # stiffness matrix  L = Inc*W*Inc'
        # mass matrix       M = 1/6*(abs.(Inc*W_inv*Inc')+abs.(Diagonal(Inc*W_inv*Inc')))
        # return L,M
        return SparseMatrixCSC(Inc*W*Inc'),SparseMatrixCSC(1/6*(abs.(Inc*W_inv*Inc')+abs.(Diagonal(Inc*W_inv*Inc'))))
    else 
        # stiffness matrix  L = Inc*W*Inc'
        # mass matrix       M = Diagonal(1/2*diag(Inc*W_inv*Inc'))
        # return L,M 
        return SparseMatrixCSC(Inc*W*Inc'),Diagonal(1/2*diag(Inc*W_inv*Inc'))
    end
end


""" left hat function """
lefthat(x,k,h)=1-(k*h-x)/h
""" right hat function """
righthat(x,k,h)=1+(k*h-x)/h
""" vertex hat function if v is start vertex """
vertexhatstart(x,h)=1-(x/h)
""" vertex hat function if v is end vertex """
vertexhatend(x,h,ℓ)=1-(ℓ-x)/h

"""
    fe_rhs(Γ::MetricGraph, rhs::Vector{Function}, h_max::Number)

Return discretized right-hand side rhs with maximum step size h_max on each edge. 
"""
function fe_rhs(Γ::MetricGraph, rhs::Vector{Function}, h_max::Number)
        u_V = zeros(nv(Γ.G)); u_E = []
        for (j,e) in enumerate(edges(Γ.G))
            Nx_e = Int(ceil(Γ.ℓ_vec[j]/h_max)); h = Γ.ℓ_vec[j]/Nx_e;
            u_e = zeros(Int(Nx_e-1))
            ## vertex hat function start
            objectfuncstart(x) = vertexhatstart(x,h)*rhs[j](x)
            u_V[src(e)]+=quadgk(objectfuncstart,0,h)[1]
            # inner hat functions
            for k=1:Nx_e-1
                # left hat
                objectfuncleft(x) = lefthat(x,k,h)*rhs[j](x)
                u_e[k]+=quadgk(objectfuncleft,(k-1)*h,k*h)[1]
                # right hat
                objectfuncright(x) = righthat(x,k,h)*rhs[j](x)
                u_e[k]+=quadgk(objectfuncright,k*h,(k+1)*h)[1]
            end
            # vertex hat function end
            objectfuncend(x) = vertexhatend(x,h,Nx_e*h)*rhs[j](x)
            u_V[dst(e)]+=quadgk(objectfuncend,Nx_e*h-h,Nx_e*h)[1]
            u_E=vcat(u_E,u_e)
        end
        return vcat(u_V,u_E)
    end

"""
    fe_rhs(Γ::EquilateralMetricGraph, rhs::Vector{Function}, h_max::Number)

Simplified version for equilateral graphs.
"""
function fe_rhs(Γ::EquilateralMetricGraph, rhs::Vector{Function}, h_max::Number)
    Nx = Int(ceil(Γ.ℓ/h_max)); h = Γ.ℓ/Nx;
    u_V = zeros(nv(Γ.G)); u_E = zeros(Int.(ne(Γ.G)*Nx-ne(Γ.G)))
    e_count = 0
    for (j,e) in enumerate(edges(Γ.G))
        ## vertex hat function start
        objectfuncstart(x) = vertexhatstart(x,h)*rhs[j](x)
        u_V[src(e)]+=quadgk(objectfuncstart,0,h)[1]
        # inner hat functions
        for k=1:Nx-1
            # left hat
            objectfuncleft(x) = lefthat(x,k,h)*rhs[j](x)
            u_E[e_count+k]+=quadgk(objectfuncleft,(k-1)*h,k*h)[1]
            # right hat
            objectfuncright(x) = righthat(x,k,h)*rhs[j](x)
            u_E[e_count+k]+=quadgk(objectfuncright,k*h,(k+1)*h)[1]
            # together
        end
        # vertex hat function end
        objectfuncend(x) = vertexhatend(x,h,Nx*h)*rhs[j](x)
        u_V[dst(e)]+=quadgk(objectfuncend,Nx*h-h,Nx*h)[1]
        e_count+=Nx-1
    end
    return vcat(u_V,u_E)
end

struct FiniteElementDiscretizationElliptic
    # stiffness matrix
    stiff::SparseMatrixCSC
    # mass matrix
    mass::SparseMatrixCSC
    # right hand side
    f::Vector
    # inner grid points per edge
    Nx::Union{Int,Vector} 
end

struct FiniteElementDiscretizationHeat
    # stiffness matrix
    stiff::SparseMatrixCSC
    # mass matrix
    mass::SparseMatrixCSC
    # inner grid points per edge
    Nx::Union{Int,Vector} 
end

"""
    fe_discretization(TP::EllipticTestProblem, J::Int; mass_aprox = false)

Compute finite element discretization of elliptic test problem 'TP' with step size 'J'. 
"""
function fe_discretization(TP::EllipticTestProblem, J::Int; mass_aprox = false)
    h_max = 2.0^(-J)
    L,M = fe_matrices(TP.Γ, h_max, mass_approx = mass_aprox);
    f   = fe_rhs(TP.Γ, TP.rhs, h_max);
    Nx  = fe_discretization_parameter(TP.Γ, h_max)
    return FiniteElementDiscretizationElliptic(L, M, f, Nx)
end

"""
    fe_discretization(TP::TPGeneralizedHeat, J::Int; mass_aprox = false)

Compute finite element discretization of generalized heat equation 'TP' with step size 2^(-'J'). 
"""
function fe_discretization(TP::TPGeneralizedHeat, J::Int; mass_aprox = false)
    h_max = 2.0^(-J)
    L,M = fe_matrices(TP.Γ, h_max, mass_approx = mass_aprox);
    Nx  = fe_discretization_parameter(TP.Γ, h_max)
    return FiniteElementDiscretizationHeat(L, M, Nx)
end

"""
    finite_element_solver(TP::EllipticTestProblem, J::Int; solver = "backslash")

Solve elliptic test problem 'TP' using a finite element discretiation with maximum step size 2^(-'J').

The backslach operator is used as a default solver for the semidiscretized system. Set solver = "multigrid" to apply multigrid solver.
"""
function finite_element_solver(TP::EllipticTestProblem, J::Int; solver = "backslash")
    if solver == "backslah"
        h_max = 2.0^(-J)
        FE = fe_discretization(TP, h_max);
        return u_fe = (FE.stiff+TP.pot*FE.mass)\FE.f;
    elseif solver == "multigrid"
        return ni_mgm(TP, J)
    end
end

"""
    finite_element_solver(TP::TPGeneralizedHeat, T::Number, J::Int; solver = "multigrid")

Solve generalized heat equation test problem 'TP' at time 'T' using a finite element discretiation with maximum step size 2^(-'J') followed by CN-MGM.

"""
function finite_element_solver(TP::TPGeneralizedHeat, T::Number, J::Int; solver = "multigrid")
    if solver == "multigrid"
        return cn_mgm(TP, T, J)
    end
end

"""
    fe_ltwo_error(Γ::EquilateralMetricGraph, u_fe::Vector, u::Vector{Function})

Compute L_2 error of the finite element solution with coefficients u_fe with 
respect to exact solution u
"""
function fe_ltwo_error(Γ::EquilateralMetricGraph, u_fe::Vector, u::Vector{Function})
    n = nv(Γ.G); m = ne(Γ.G); NxE = Int((length(u_fe)-n)/m+1); h = Γ.ℓ/NxE;
    # initialize l2_error
    l2_error_sqr = 0; 
    # counters for vector access
    e_count = n; 
    for (m_count,e) in enumerate(edges(Γ.G))
        # assemble coefficient vector for current edge
        uh_e = zeros(NxE+1)
        uh_e[1] = u_fe[src(e)]; uh_e[NxE+1] = u_fe[dst(e)]
        uh_e[2:end-1] = u_fe[e_count+1:e_count+NxE-1]
        e_count+=NxE-1
        # assemble spline
        alloc_points_h = 0:h:ℓ
        basis_h = BSplineBasis(2,alloc_points_h)
        spl_h = Spline(basis_h,uh_e)
        ## L2-error
        for k in 0:NxE-1
            Ltwointegrand(x) = ((splinevalue(spl_h,x))-u_star[m_count](x))^2
            l2_error_sqr+=quadgk(Ltwointegrand,k*h,(k+1)*h)[1]
        end
        
    end
    # return norm
    return sqrt(l2_error_sqr)
end


"""
    fe_error(Γ::MetricGraph, u_fe::Vector, u::Vector{Function}, u_deriv::Vector{Function}, Nx_vec::Vector)

Compute L_2 and H^1 error of the finite element solution with coefficients u_fe with 
respect to exact solution u with derivative u_deriv.
"""
function fe_error(Γ::MetricGraph, u_fe::Vector, u::Vector{Function}, u_deriv::Vector{Function}, Nx_vec::Vector)
    n=nv(Γ.G); m=ne(Γ.G); 
    
    l2_error_sqr=0; h1_error_sqr=0
    e_count=n
    for (j,e) in enumerate(edges(Γ.G))
        NxE = Nx_vec[j]
        h = Γ.ℓ_vec[j]/NxE
        ## current solution
        # assemble coefficient vector for current edge
        uh_e=zeros(NxE+1)
        uh_e[1]=u_fe[src(e)]; uh_e[NxE+1]=u_fe[dst(e)]
        uh_e[2:end-1]=u_fe[e_count+1:e_count+NxE-1]
        e_count+=NxE-1
        # assemble spline
        alloc_points_h=0:h:Γ.ℓ_vec[j]
        basis_h=BSplineBasis(2,alloc_points_h)
        spl_h = Spline(basis_h,uh_e)

        ## L2-error
        for k in 0:NxE-1
            Ltwointegrand(x)=(splinevalue(spl_h,x)-u[j](x))^2
            l2_error_sqr+=quadgk(Ltwointegrand,k*h,(k+1)*h)[1]
        end
        Honeintegrand(x)=(splinevalue(spl_h,x,Derivative(1))-u_deriv[j](x))^2
        # H1-error
        for k in 0:NxE-1
            h1_error_sqr+=quadgk(Honeintegrand, k*h,(k+1)*h)[1]
        end

    end
    return [sqrt(l2_error_sqr),sqrt(l2_error_sqr+h1_error_sqr)]
end

"""
    fe_error(Γ::EquilateralMetricGraph, u_fe::Vector, u::Vector{Function}, u_deriv::Vector{Function})

Equilateral version.
"""
function fe_error(Γ::EquilateralMetricGraph, u_fe::Vector, u::Vector{Function}, u_deriv::Vector{Function})
    n=nv(Γ.G); m=ne(Γ.G); NxE=Int((length(u_fe)-n)/m+1); h=Γ.ℓ/NxE; 
    l2_error_sqr=0; h1_error_sqr=0
    e_count=n
    for (j,e) in enumerate(edges(Γ.G))

        ## current solution
        # assemble coefficient vector for current edge
        uh_e=zeros(NxE+1)
        uh_e[1]=u_fe[src(e)]; uh_e[NxE+1]=u_fe[dst(e)]
        uh_e[2:end-1]=u_fe[e_count+1:e_count+NxE-1]
        e_count+=NxE-1
        # assemble spline
        alloc_points_h=0:h:Γ.ℓ
        basis_h=BSplineBasis(2,alloc_points_h)
        spl_h = Spline(basis_h,uh_e)

        ## L2-error
        for k in 0:NxE-1
            Ltwointegrand(x)=(splinevalue(spl_h,x)-u[j](x))^2
            l2_error_sqr+=quadgk(Ltwointegrand,k*h,(k+1)*h)[1]
        end
        Honeintegrand(x)=(splinevalue(spl_h,x,Derivative(1))-u_deriv[j](x))^2
        # H1-error
        for k in 0:NxE-1
            h1_error_sqr+=quadgk(Honeintegrand, k*h,(k+1)*h)[1]
        end

    end
    return [sqrt(l2_error_sqr),sqrt(l2_error_sqr+h1_error_sqr)]
end


"""
    fe_error(TP::EllipticTestProblem, u_fe::Vector)

Compute L_2 and H^1 error of the finite element solution with coefficients 'u_fe' with 
respect to exact solution 'TP.u' with derivative 'TP.u_deriv'.
"""
function fe_error(TP::EllipticTestProblem, u_fe::Vector, Nx::Union{Int,Vector})
    n = nv(TP.Γ.G); m = ne(TP.Γ.G); 
    l2_error_sqr = 0; h1_error_sqr = 0;
    e_count = n;
    if typeof(TP.Γ) == EquilateralMetricGraph
        Nx = Nx*ones(Int,m)
        ℓ_vec = TP.Γ.ℓ * ones(m)
    else
        ℓ_vec = TP.Γ.ℓ_vec
    end
    for (j,e) in enumerate(edges(TP.Γ.G))
        Nxe = Nx[j]; h = ℓ_vec[j]/Nxe
        ## current solution
        # assemble coefficient vector for current edge
        uh_e = zeros(Nxe+1)
        uh_e[1] = u_fe[src(e)]; uh_e[Nxe+1] = u_fe[dst(e)]
        uh_e[2:end-1] = u_fe[e_count+1:e_count+Nxe-1]
        e_count += Nxe-1
        # assemble spline
        alloc_points_h = 0:h:ℓ_vec[j]
        basis_h = BSplineBasis(2,alloc_points_h)
        spl_h = Spline(basis_h,uh_e)

        ## L2-error
        for k in 0:Nxe-1
            Ltwointegrand(x) = (splinevalue(spl_h,x)-TP.u[j](x))^2
            l2_error_sqr += quadgk(Ltwointegrand,k*h,(k+1)*h)[1]
        end
        Honeintegrand(x) = (splinevalue(spl_h,x,Derivative(1))-TP.u_deriv[j](x))^2
        # H1-error
        for k in 0:Nxe-1
            h1_error_sqr += quadgk(Honeintegrand, k*h,(k+1)*h)[1]
        end

    end
    return [sqrt(l2_error_sqr),sqrt(l2_error_sqr+h1_error_sqr)]
end
