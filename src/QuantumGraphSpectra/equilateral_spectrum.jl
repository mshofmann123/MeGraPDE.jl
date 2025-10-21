# © Anna Weller, University of Cologne, 2023

"""
    harmonic_eigen(G::SimpleGraph)

Return eigenvalues and eigenvectors of the harmonic graph Laplacian Matrix ``\\Delta_{\\mathcal{G}}``

A complete eigendecomposition of the normalized graph Laplacian is computed using julias "eigen". 
The harmonic graph Laplacian is similar to the normalized Laplacian, in particular the eigenvalues 
``\\mu`` agree and ``\\Phi = \\mathbf{D}^{-\\frac{1}{2}} \\Psi``
"""
function harmonic_eigen(G::SimpleGraph)
    # Laplacian and Degree matrix
    L = laplacian_matrix(G)
    D_sqrt_inv = Diagonal(diag(L) .^ (-1/2))
    # eigendecomposition of L_norm 
    Mu_vec, Ψ_vec=eigen(Symmetric(Matrix(D_sqrt_inv*L*D_sqrt_inv)))
    # return vector containing eigenvalues and matrix containing eigenvectors [EV_1,...,EV_n] 
    return Mu_vec, D_sqrt_inv*Ψ_vec # D^{-1/2}*Ψ_vec are the eigenvectors of L_harm 
end

"""
    eigvals_quantum(Γ::EquilateralMetricGraph; K=3, sorted=true, only_vertex=false)

Compute all eigenvalues ``\\lambda < ('K' \\pi)/ \\ell)^2`` of the equilateral metric graph 'Γ'.
"""
function eigvals_quantum(Γ::EquilateralMetricGraph; K=3, sorted=true, only_vertex=false)
    # eigvals of harmonic graph Laplacian
    μ_vec, _ = harmonic_eigen(Γ.G)
    # determine number of vertex and non-vertex eigvals
    n=nv(Γ.G);
    m=ne(Γ.G)
    if is_bipartite(Γ.G)
        no_vertex = n-2;
        no_non_vertex_even = m-n+2;
        no_non_vertex_odd = m-n+2
    else
        no_vertex = n-1;
        no_non_vertex_even = m-n+2;
        no_non_vertex_odd = m-n
    end

    Λ = []
    for k in 0:(K - 1)
        # k even
        if k%2==0
            ## non-vertex eigenvalues, k even
            if only_vertex == false
                if k==0
                    push!(Λ, 0)
                else
                    for _ in 1:no_non_vertex_even
                        push!(Λ, (k*pi/Γ.ℓ)^2)
                    end
                end
            end
            # vertex eigenvalues, k even
            for μ in μ_vec[2:(no_vertex + 1)]
                λ = ((1/Γ.ℓ)*(acos(1-μ)+k*pi))^2
                push!(Λ, λ)
            end
            # k odd
        else
            ## non-vertex eigenvalues, k odd
            if only_vertex == false
                for _ in 1:no_non_vertex_odd
                    push!(Λ, (k*pi/Γ.ℓ)^2)
                end
            end
            ## vertex eigevanlues, k odd
            if sorted == true
                for μ in reverse(μ_vec[2:(no_vertex + 1)])
                    λ = ((1/Γ.ℓ)*(-acos(1-μ)+(k+1)*pi))^2
                    push!(Λ, λ)
                end
            else
                for μ in (μ_vec[2:(no_vertex + 1)])
                    λ = ((1/Γ.ℓ)*(-acos(1-μ)+(k+1)*pi))^2
                    push!(Λ, λ)
                end
            end
        end
    end
    # return
    if only_vertex==true
        return QuantumGraphVertexEigvals(length(Λ), Λ)
    else
        return QuantumGraphEigvals(length(Λ), Λ)
    end
end

"""
    norm_equilateral_eigenfunction(ϕ::EquilateralEigfunc)

Compute norm of eigenfunction 'ϕ' of equilateral metric graph
"""
function eigfunc_norm(λ::Number, A::Vector, B::Vector, ℓ::Number, m::Int)
    norm = 0
    λ_sqrt = sqrt(λ)
    for j in 1:m
        Ae = A[j];
        Be = B[j]
        norm +=
            (
                2*λ_sqrt*ℓ*(Ae^2+Be^2)+(Ae^2-Be^2)*sin(2*λ_sqrt*ℓ)+2*Ae*Be*(
                    1-cos(2*λ_sqrt*ℓ)
                )
            )/(4*λ_sqrt)
    end
    return sqrt(norm)
end

"""
    norm_equilateral_eigenfunction(ϕ::EquilateralEigfunc)

Compute norm of eigenfunction 'ϕ' of equilateral metric graph
"""
function eigfunc_norm_vertex(Γ::EquilateralMetricGraph, λ::Number, μ::Number, Φ::Vector)
    norm = 0
    λ_sqrt = sqrt(λ)
    for (j, e) in enumerate(edges(Γ.G))
        v1 = src(e);
        v2 = dst(e)
        norm += 1/2*(Φ[v1]^2 + (1/sin(acos(1-μ))*(Φ[v2]-Φ[v1]*(1-μ)))^2)*Γ.ℓ
        norm +=
            1/(4*λ_sqrt) *
            (Φ[v1]^2 - (1/sin(acos(1-μ))*(Φ[v2]-Φ[v1]*(1-μ)))^2) *
            sin(2*acos(1-μ))
        norm +=
            1/(4*λ_sqrt) *
            (2*Φ[v1]*(1/sin(acos(1-μ))*(Φ[v2]-Φ[v1]*(1-μ)))*(1-cos(2*acos(1-μ))))
    end
    return sqrt(norm)
end

"""
    eigen_quantum(Γ::EquilateralMetricGraph; K=3, sorted=true, sparse_svd=false)

Compute all eigenvalues ``\\lambda < ('K' \\pi)/ \\ell)^2`` and corresponding eigenfunctions ``\\phi`` with 
```math
{\\phi}_e = A_e \\cos(\\sqrt{\\lambda} x) + B_e \\sin (\\sqrt{\\lambda} x)
```
of the equilateral metric graph 'Γ'.

The coefficient ``A_e``, ``B_e`` are stored in A = [Ae1,…,Aem]' and B = [Be1,…,Bem]' for 
each eigenfunction. The coefficients are normalized such that all eigenfunctions fulfill `` \\| \\phi \\|=1 ``.
"""
function eigen_quantum(Γ::EquilateralMetricGraph; K=3, sparse_svd=false, round_off=false)
    n = nv(Γ.G);
    m = ne(Γ.G)

    # compute harmonic eigendecompositon
    Mu_vec, Φ_vec = harmonic_eigen(Γ.G)
    if round_off == true
        Φ_vec = map!(x -> isapprox(x, 0; atol=10^(-15)) ? 0 : x, Φ_vec, Φ_vec)
    end
    # determine number of vertex and non-vertex eigenvalues
    if is_bipartite(Γ.G)
        no_vertex_ev = n-2
        no_non_vertex_ev_odd = m-n+2
        no_non_vertex_ev_even = m-n+2
    else
        no_vertex_ev = n-1
        no_non_vertex_ev_odd = m-n
        no_non_vertex_ev_even = m-n+2
    end

    # determine total number of eigenvalues for K
    if K%2 == 0
        total_no_evals = Int(
            1+K*no_vertex_ev+(K/2)*no_non_vertex_ev_odd+(K/2-1)*no_non_vertex_ev_even
        )
    else
        total_no_evals = Int(
            1+K*no_vertex_ev+((K-1)/2)*no_non_vertex_ev_odd+((K-1)/2)*no_non_vertex_ev_even
        )
    end

    # assemble output matrices
    Λ = zeros(total_no_evals);
    A = spzeros(m, total_no_evals);
    B = spzeros(m, total_no_evals)

    eval_count = 1
    for k in 0:(K - 1)
        # k even 
        if k%2 == 0
            if k == 0
                norm = (sqrt(m*Γ.ℓ))
                Λ[eval_count] = 0;
                A[:, eval_count] = ones(m)/norm
                # next eval
                eval_count+=1
            else
                ### non-vertex k even ##################################
                L̃ = extended_laplacian(Γ, k)
                sqrtD⁻¹ = Diagonal(diag(L̃) .^ (-1/2))
                # quantum graph eigenvalue
                λ = ((k*pi)/Γ.ℓ)^2
                # discrete eigenvalue of extended graph
                μ = 1-cos((k/(k+1))*pi)
                # compute svd
                if sparse_svd == false
                    _, _, V = svd(Matrix(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I))
                else
                    _, V = eigs(
                        Symmetric(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I);
                        sigma=10^(-7),
                        nev=no_non_vertex_ev_even,
                    )
                end
                for i in 1:no_non_vertex_ev_even
                    # extended eigenvectors
                    Φ_ext = sqrtD⁻¹*(V[:, end - i + 1])
                    # store coefficients
                    A_q = zeros(m);
                    B_q = zeros(m)
                    for (j, e) in enumerate(edges(Γ.G))
                        e_0 = src(e);
                        e_ℓ = dst(e)
                        w1 = n+k*(j-1)+1
                        A_q[j] = Φ_ext[e_0]
                        B_q[j] =
                            (1/sin((k/(k+1))*pi))*(Φ_ext[w1]-Φ_ext[e_0]*cos((k/(k+1))*pi))
                    end
                    # norm
                    norm = eigfunc_norm(λ, A_q, B_q, Γ.ℓ, m)
                    # store normed coefficients
                    Λ[eval_count] = λ;
                    A[:, eval_count] = A_q ./ norm;
                    B[:, eval_count] = B_q ./ norm
                    # next eval
                    eval_count += 1
                end
            end
            ### vertex k even ######################################
            for i in 2:(no_vertex_ev + 1)
                # quantum graph eigenvalue
                λ = (1/Γ.ℓ*(acos(1-Mu_vec[i])+k*pi))^2;
                Φ = Φ_vec[:, i]
                # store coefficients
                A_q = zeros(m);
                B_q = zeros(m)
                for (j, e) in enumerate(edges(Γ.G))
                    e_0 = src(e);
                    e_ℓ = dst(e)
                    A_q[j] = Φ[e_0]
                    B_q[j] = 1/(sin(acos(1-Mu_vec[i])))*(Φ[e_ℓ]-Φ[e_0]*(1-Mu_vec[i]))
                end
                # norm
                norm=eigfunc_norm_vertex(Γ, λ, Mu_vec[i], Φ)
                # store normed coefficients
                Λ[eval_count] = λ;
                A[:, eval_count] = A_q ./ norm;
                B[:, eval_count] = B_q ./ norm
                # next eval
                eval_count+=1
            end
            ### k odd ##############################################################################
        else
            ### non-vertex k odd ##################################
            L̃ = extended_laplacian(Γ, k)
            sqrtD⁻¹ = Diagonal(diag(L̃) .^ (-1/2))
            # quantum graph eigenvalue
            λ=((k*pi)/Γ.ℓ)^2
            # discrete eigenvalue of extended graph
            μ=1-cos((k/(k+1))*pi)
            # compute svd
            if sparse_svd == false
                _, _, V = svd(Matrix(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I))
            else
                _, V = eigs(
                    Symmetric(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I);
                    sigma=10^(-7),
                    nev=no_non_vertex_ev_even,
                )
            end
            #
            for i in 1:no_non_vertex_ev_odd
                Φ_ext=sqrtD⁻¹*(V[:, end - i + 1])
                # store coefficients
                A_q=zeros(m);
                B_q=zeros(m)
                for (j, e) in enumerate(edges(Γ.G))
                    e_0=src(e);
                    e_ℓ=dst(e)
                    w1=n+k*(j-1)+1
                    A_q[j]=Φ_ext[e_0]
                    B_q[j]=(1/sin((k/(k+1))*pi))*(Φ_ext[w1]-Φ_ext[e_0]*cos((k/(k+1))*pi))
                end
                # norm
                norm=eigfunc_norm(λ, A_q, B_q, Γ.ℓ, m)
                # store normed coefficients
                Λ[eval_count]=λ;
                A[:, eval_count]=A_q ./ norm;
                B[:, eval_count]=B_q ./ norm
                # next eval
                eval_count+=1
            end
            ### vertex k odd #######################################
            for i in (no_vertex_ev + 1):-1:2
                # quantum graph eigenvalue
                λ=(1/Γ.ℓ*(-acos(1-Mu_vec[i])+(k+1)*pi))^2;
                Φ=Φ_vec[:, i]
                # store coefficients
                A_q=zeros(m);
                B_q=zeros(m)
                for (j, e) in enumerate(edges(Γ.G))
                    e_0=src(e);
                    e_ℓ=dst(e)
                    A_q[j]=Φ[e_0]
                    B_q[j]=1/(sin(-acos(1-Mu_vec[i])))*(Φ[e_ℓ]-Φ[e_0]*(1-Mu_vec[i]))
                end
                #norm=eigfunc_norm(λ,A_q,B_q,Γ.ℓ,m)
                norm=eigfunc_norm_vertex(Γ, λ, Mu_vec[i], Φ)
                # store normed coefficients
                Λ[eval_count]=λ;
                A[:, eval_count]=A_q ./ norm;
                B[:, eval_count]=B_q ./ norm
                # next eval
                eval_count+=1
            end
        end
    end
    #if round_off == true
    #    B = map!(x -> isapprox(x, 0, atol=10^(-15)) ? 0 : x, B, B)
    #end
    return QuantumGraphEigen(length(Λ), Λ, sparse(A), sparse(B))
end

"""
    eigen_quantum_old(Γ::EquilateralMetricGraph; K=3, sorted=true, sparse_svd=false)

Compute all eigenvalues λ < ('K'*π/ℓ)^2 and corresponding eigenfunctions ϕ with 
ϕ_e = A_e cos(sqrt(λ)x) + B_e sin (sqrt(λ)x) of the equilateral metric graph 'Γ'.

The coefficient A_e, B_e are stored in A = [A_e1,…,A_em]' and B = [B_e1,…,B_em]' for 
each eigenfunction. The coefficients are normalized such that all eigenfunctions fulfill ||ϕ||=1.
"""
function eigen_quantum_old(Γ::EquilateralMetricGraph; K=3, sorted=true, sparse_svd=false)
    n = nv(Γ.G);
    m = ne(Γ.G)

    # compute harmonic eigendecompositon
    Mu_vec, Φ_vec = harmonic_eigen(Γ.G)

    # determine number of vertex and non-vertex eigenvalues
    if is_bipartite(Γ.G)
        no_vertex_ev = n-2
        no_non_vertex_ev_odd = m-n+2
        no_non_vertex_ev_even = m-n+2
    else
        no_vertex_ev = n-1
        no_non_vertex_ev_odd = m-n
        no_non_vertex_ev_even = m-n+2
    end

    # determine total number of eigenvalues for K
    if K%2 == 0
        total_no_evals = Int(
            1+K*no_vertex_ev+(K/2)*no_non_vertex_ev_odd+(K/2-1)*no_non_vertex_ev_even
        )
    else
        total_no_evals = Int(
            1+K*no_vertex_ev+((K-1)/2)*no_non_vertex_ev_odd+((K-1)/2)*no_non_vertex_ev_even
        )
    end

    # assemble output matrices
    Λ = zeros(total_no_evals);
    A = spzeros(m, total_no_evals);
    B = spzeros(m, total_no_evals)

    eval_count = 1
    for k in 0:(K - 1)
        # k even 
        if k%2 == 0
            if k == 0
                norm = (sqrt(m*Γ.ℓ))
                Λ[eval_count] = 0;
                A[:, eval_count] = ones(m)/norm
                # next eval
                eval_count+=1
            else
                ### non-vertex k even ##################################
                L̃ = extended_laplacian(Γ, k)
                sqrtD⁻¹ = Diagonal(diag(L̃) .^ (-1/2))
                # quantum graph eigenvalue
                λ = ((k*pi)/Γ.ℓ)^2
                # discrete eigenvalue of extended graph
                μ = 1-cos((k/(k+1))*pi)
                # compute svd
                if sparse_svd == false
                    _, _, V = svd(Matrix(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I))
                else
                    _, V = eigs(
                        Symmetric(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I);
                        sigma=10^(-7),
                        nev=no_non_vertex_ev_even,
                    )
                end
                for i in 1:no_non_vertex_ev_even
                    # extended eigenvectors
                    Φ_ext = sqrtD⁻¹*(V[:, end - i + 1])
                    # store coefficients
                    A_q = zeros(m);
                    B_q = zeros(m)
                    for (j, e) in enumerate(edges(Γ.G))
                        e_0 = src(e);
                        e_ℓ = dst(e)
                        w1 = n+k*(j-1)+1
                        A_q[j] = Φ_ext[e_0]
                        B_q[j] =
                            (1/sin((k/(k+1))*pi))*(Φ_ext[w1]-Φ_ext[e_0]*cos((k/(k+1))*pi))
                    end
                    # norm
                    norm = eigfunc_norm(λ, A_q, B_q, Γ.ℓ, m)
                    # store normed coefficients
                    Λ[eval_count] = λ;
                    A[:, eval_count] = A_q ./ norm;
                    B[:, eval_count] = B_q ./ norm
                    # next eval
                    eval_count += 1
                end
            end
            ### vertex k even ######################################
            for i in 2:(no_vertex_ev + 1)
                # quantum graph eigenvalue
                λ = (1/Γ.ℓ*(acos(1-Mu_vec[i])+k*pi))^2;
                Φ = Φ_vec[:, i]
                # store coefficients
                A_q = zeros(m);
                B_q = zeros(m)
                for (j, e) in enumerate(edges(Γ.G))
                    e_0 = src(e);
                    e_ℓ = dst(e)
                    A_q[j] = Φ[e_0]
                    B_q[j] = 1/(sin(sqrt(λ)*Γ.ℓ))*(Φ[e_ℓ]-Φ[e_0]*cos(sqrt(λ)*Γ.ℓ))
                end
                # norm
                norm=eigfunc_norm(λ, A_q, B_q, Γ.ℓ, m)
                # store normed coefficients
                Λ[eval_count] = λ;
                A[:, eval_count] = A_q ./ norm;
                B[:, eval_count] = B_q ./ norm
                # next eval
                eval_count+=1
            end
            ### k odd ##############################################################################
        else
            ### non-vertex k odd ##################################
            L̃ = extended_laplacian(Γ, k)
            sqrtD⁻¹ = Diagonal(diag(L̃) .^ (-1/2))
            # quantum graph eigenvalue
            λ=((k*pi)/Γ.ℓ)^2
            # discrete eigenvalue of extended graph
            μ=1-cos((k/(k+1))*pi)
            # compute svd
            if sparse_svd == false
                _, _, V = svd(Matrix(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I))
            else
                _, V = eigs(
                    Symmetric(sqrtD⁻¹*L̃*sqrtD⁻¹-μ*I);
                    sigma=10^(-7),
                    nev=no_non_vertex_ev_even,
                )
            end
            #
            for i in 1:no_non_vertex_ev_odd
                Φ_ext=sqrtD⁻¹*(V[:, end - i + 1])
                # store coefficients
                A_q=zeros(m);
                B_q=zeros(m)
                for (j, e) in enumerate(edges(Γ.G))
                    e_0=src(e);
                    e_ℓ=dst(e)
                    w1=n+k*(j-1)+1
                    A_q[j]=Φ_ext[e_0]
                    B_q[j]=(1/sin((k/(k+1))*pi))*(Φ_ext[w1]-Φ_ext[e_0]*cos((k/(k+1))*pi))
                end
                # norm
                norm=eigfunc_norm(λ, A_q, B_q, Γ.ℓ, m)
                # store normed coefficients
                Λ[eval_count]=λ;
                A[:, eval_count]=A_q ./ norm;
                B[:, eval_count]=B_q ./ norm
                # next eval
                eval_count+=1
            end
            ### vertex k odd #######################################
            for i in (no_vertex_ev + 1):-1:2
                # quantum graph eigenvalue
                λ=(1/Γ.ℓ*(acos(1-Mu_vec[i])-(k+1)*pi))^2;
                Φ=Φ_vec[:, i]
                # store coefficients
                A_q=zeros(m);
                B_q=zeros(m)
                for (j, e) in enumerate(edges(Γ.G))
                    e_0=src(e);
                    e_ℓ=dst(e)
                    A_q[j]=Φ[e_0]
                    B_q[j]=1/(sin(sqrt(λ)*Γ.ℓ))*(Φ[e_ℓ]-Φ[e_0]*cos(sqrt(λ)*Γ.ℓ))
                end
                # norm
                norm=eigfunc_norm(λ, A_q, B_q, Γ.ℓ, m)
                # store normed coefficients
                Λ[eval_count]=λ;
                A[:, eval_count]=A_q ./ norm;
                B[:, eval_count]=B_q ./ norm
                # next eval
                eval_count+=1
            end
        end
    end

    return QuantumGraphEigen(length(Λ), Λ, sparse(A), sparse(B))
end

"""
    eigenfunction(Γ::EquilateralMetricGraph, σ::QuantumGraphEigen, q::Int)

Return 'q'th eigenfunction in 'σ' (closed form).
"""
function eigenfunction(Γ::EquilateralMetricGraph, σ::QuantumGraphEigen, q::Int)
    ϕ_q = Vector{Function}(undef, ne(Γ.G));
    for j in 1:ne(Γ.G)
        ϕ_q[j] = x -> σ.A[j, q]*cos(sqrt(σ.Λ[q])*x)+σ.B[j, q]*sin(sqrt(σ.Λ[q])*x)
    end
    return ϕ_q
end

"""
    count_eigvals_K(Γ::EquilateralMetricGraph, K::Int)

Return number of eigenvalues `` \\lambda < ('K'* \\pi/ \\ell)^2 ``.

"""
function count_eigvals_K(Γ::EquilateralMetricGraph, K::Int)
    # determine number of vertex and non-vertex eigenvalues
    n = nv(Γ.G);
    m = ne(Γ.G)
    if is_bipartite(Γ.G)
        no_vertex_ev = n-2
        no_non_vertex_ev_odd = m-n+2
        no_non_vertex_ev_even = m-n+2
    else
        no_vertex_ev = n-1
        no_non_vertex_ev_odd = m-n
        no_non_vertex_ev_even = m-n+2
    end
    # k_even
    if K%2 == 0
        return Int(
            1 + K*no_vertex_ev + K/2*no_non_vertex_ev_odd + (K/2-1)*no_non_vertex_ev_even
        )
    else
        return Int(
            1 +
            K*no_vertex_ev +
            (K-1)/2*no_non_vertex_ev_even +
            ((K+1)/2-1)*no_non_vertex_ev_odd,
        )
    end
end
