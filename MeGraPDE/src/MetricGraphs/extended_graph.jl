# © Anna Weller, University of Cologne, 2023

"""
    extended_incidence_matrix(Γ::MetricGraph, h_max::Number)

Return extended incidence matrix of 'Γ' with maximal step length 'h_max' per edge.

Construction of incidence matrix via kron-products according to (AB) (section 4.1), see also (W), section 3.1 for a summary.
"""
function extended_incidence_matrix(Γ::MetricGraph, h_max::Number)
    # store step lengths on each edge
    h_vec = zeros(length(Γ.ℓ_vec))
    # assemble in and out incidence matrix of G
    bfN = incidence_matrix(Γ.G,oriented=true)
    N_in = 0.5*(bfN+abs.(bfN)); N_out = 0.5*(bfN-abs.(bfN))
    # assemble N_tilde 
    N_tilde_V = spzeros(nv(Γ.G),0)
    N_tilde_E = spzeros(0,0)
    # iterate over edges
    for j = 1:ne(Γ.G) 
        # number of grid points on current edge
        Nxe = Int.(ceil(Γ.ℓ_vec[j]/h_max)); h_vec[j] = Γ.ℓ_vec[j]/Nxe
        # comput (N_tilde_V)_e
        bfe1 = I(Nxe)[1,:]; bfeNxe = I(Nxe)[Nxe,:]
        N_tilde_Ve = kron(N_out[:,j],bfe1')+kron(N_in[:,j],bfeNxe')
        N_tilde_V = hcat(N_tilde_V,N_tilde_Ve)
        # compute (N_tilde_E)_e
        N_tilde_Ee = spdiagm(Nxe-1,Nxe,0 => ones(Nxe-1),1 => -ones(Nxe-1))
        N_tilde_E = blockdiag(N_tilde_E,N_tilde_Ee)
    end
    # return
    return SparseMatrixCSC(vcat(N_tilde_V,N_tilde_E))
end

"""
    extended_incidence_matrix(Γ::EquilateralMetricGraph, h_max::Number)

Equilateral version with some simplifications.
"""
function extended_incidence_matrix(Γ::EquilateralMetricGraph, h_max::Number)
    # number of grid points on each edge
    Nx_vec = Int(ceil(Γ.ℓ/h_max)); 
    # assemble in and out incidence matrix of G
    bfN = incidence_matrix(Γ.G,oriented=true)
    N_in = 0.5*(bfN+abs.(bfN)); N_out = 0.5*(bfN-abs.(bfN))
    # assemble N_tilde 
    N_tilde_V = spzeros(nv(Γ.G),0)
    N_tilde_E = spzeros(0,0)
    # incidence matrix of edges
    N_tilde_Ee = spdiagm(Nx_vec-1, Nx_vec,0 => ones(Nx_vec-1), 1 => -ones(Nx_vec-1))
    # iterate over edges
    for j = 1:ne(Γ.G) 
        # comput (N_tilde_V)_e
        bfe1 = I(Nx_vec)[1,:]; bfeNx_vec = I(Nx_vec)[Nx_vec,:]
        N_tilde_Ve = kron(N_out[:,j],bfe1')+kron(N_in[:,j],bfeNx_vec')
        N_tilde_V = hcat(N_tilde_V,N_tilde_Ve)
        # compute (N_tilde_E)_e
        N_tilde_E = blockdiag(N_tilde_E,N_tilde_Ee)
    end
    return SparseMatrixCSC(vcat(N_tilde_V,N_tilde_E))
end

"""
    extended_laplacian(Γ::EquilateralMetricGraph, k::Int)

Compute extended graph Laplacian matrix of 'Γ' with 'k' artificial vertices on each edge.

k inner vertices means each edge is partitioned in k+1 subdivision.
The construction of the graph Laplacian relies on the same manipulations
of the original graph as in the extended_incidence_matrix routine. Here,
however, the Laplacian matrix L=NN^T is returned and some simplifications due 
to the equilateral edge length apply.
"""
function extended_laplacian(Γ::EquilateralMetricGraph, k::Int)
    # assemble in and out incidence matrix of G
    bfN = incidence_matrix(Γ.G,oriented=true)
    N_in = 0.5*(bfN+abs.(bfN)); N_out = 0.5*(bfN-abs.(bfN))
    # assemble N_tilde 
    N_tilde_V = spzeros(nv(Γ.G),0)
    N_tilde_E = spzeros(0,0)
    # for each edge we obtain incidence matrix of the form
    N_tilde_Ee = spdiagm(k,k+1,0=>ones(k),1=>-ones(k))
    for j = 1:ne(Γ.G) 
        # comput (N_tilde_V)N_tilde_E
        bfe1 = I(k+1)[1,:]; bfeNxe = I(k+1)[k+1,:]
        N_tilde_Ve = kron(N_out[:,j],bfe1')+kron(N_in[:,j],bfeNxe')
        N_tilde_V = hcat(N_tilde_V,N_tilde_Ve)
        # compute (N_tilde_E)_e
        N_tilde_E = blockdiag(N_tilde_E,N_tilde_Ee)
    end
    N = vcat(N_tilde_V,N_tilde_E)
    return SparseMatrixCSC(N*N') # graph Laplacian matrix
end


"""
    discretize_function(Γ::MetricGraph, u::Vector{Function}, h_max::Number)

Return discretized version of 'u' on the extended graph of 'Γ' with step size 'h_max' on the edges.

"""
function discretize_function(Γ::MetricGraph, u::Vector{Function}, h_max::Number)
    u_V = zeros(nv(Γ.G))
    Nx_e_sum = 0
    for (j,e) in enumerate(edges(Γ.G))
        Nx_e = Int(ceil(Γ.ℓ_vec[j]/h_max)); 
        h = Γ.ℓ_vec[j]/Nx_e
        u_V[src(e)] = u[j](0)
        u_V[dst(e)] = u[j]((Nx_e)*h)
        Nx_e_sum += Nx_e-1
    end
    u_E = zeros(Int.(Nx_e_sum))
    e_count = 0
    for j = 1:ne(Γ.G)
        Nx_e = Int(ceil(Γ.ℓ_vec[j]/h_max)); h = Γ.ℓ_vec[j]/Nx_e
        for k = 1:Nx_e-1
            u_E[e_count+k] = u[j](k*h)
        end
        e_count += Nx_e-1
    end

    return vcat(u_V,u_E)
end
    
"""
    discretize_function(Γ::MetricGraph, u::Vector{Function}, Nx_vec::Vector)

Return discretized version of 'u' on the extended graph of 'Γ' with inner grid points in 'Nx_vec'.

"""
function discretize_function(Γ::MetricGraph, u::Vector{Function}, Nx_vec::Vector)
    u_V = zeros(nv(Γ.G))
    Nx_e_sum = 0
    for (j,e) in enumerate(edges(Γ.G))
        Nx_e = Nx_vec[j]; 
        h = Γ.ℓ_vec[j]/Nx_e
        u_V[src(e)] = u[j](0)
        u_V[dst(e)] = u[j]((Nx_e)*h)
        Nx_e_sum += Nx_e-1
    end
    u_E = zeros(Int.(Nx_e_sum))
    e_count = 0
    for j = 1:ne(Γ.G)
        Nx_e = Nx_vec[j]; h = Γ.ℓ_vec[j]/Nx_e
        for k = 1:Nx_e-1
            u_E[e_count+k] = u[j](k*h)
        end
        e_count += Nx_e-1
    end

    return vcat(u_V,u_E)
end

"""
    discretize_function(Γ::EquilateralMetricGraph, u::Vector{Function}, h_max::Number)

Equilateral version.
"""
function discretize_function(Γ::EquilateralMetricGraph, u::Vector{Function}, h_max::Number)
    u_V = zeros(nv(Γ.G));
    Nx_sum = 0
    for (j,e) in enumerate(edges(Γ.G))
        Nx = ceil(Γ.ℓ/h_max); h = Γ.ℓ/Nx;
        u_V[src(e)] = u[j](0)
        u_V[dst(e)] = u[j](Nx*h)
        Nx_sum += Nx;
    end
    u_E = zeros(Int.(Nx_sum-ne(Γ.G)))
    e_count = 0
    for j = 1:ne(Γ.G)
        Nx = Int(ceil(Γ.ℓ/h_max)); h=Γ.ℓ/Nx;
        for k = 1:Nx-1
            u_E[e_count+k] = u[j](k*h)
        end
        e_count+=Nx-1
    end
    return vcat(u_V,u_E)
end