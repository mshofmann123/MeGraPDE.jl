# © Anna Weller, University of Cologne, 2023


using Graphs

"""
    spectral_solution(Γ::AbstractMetricGraph, σ::QuantumGraphEigen, coef::Vector)

Explicitly construct spectral solution u_Q on 'Γ' from eigenbasis 'σ' and coefficents 'coef'.
"""
function spectral_solution(Γ::AbstractMetricGraph, σ::QuantumGraphEigen, coef::Vector)
    Λ_sqrt = sqrt.(σ.Λ);
    u_Q = Vector{Function}(undef,ne(Γ.G));
    for j = 1:ne(Γ.G)
        function u_Q_e(x)
            Ae_Cos_x = [σ.A[j,q]*cos(Λ_sqrt[q]*x) for q=1:length(coef)];
            Be_Sin_x = [σ.B[j,q]*sin(Λ_sqrt[q]*x) for q=1:length(coef)];
            return coef'*(Ae_Cos_x+Be_Sin_x)
        end
        u_Q[j] = u_Q_e;
    end
    return u_Q
end


"""
    spectral_solution_deriv(Γ::AbstractMetricGraph, σ::QuantumGraphEigen, coef::Vector)

Explicitly construct derivative u_Q' of spectral solution on 'Γ' from eigenbasis 'σ' and coefficents 'coef'.

"""
function spectral_solution_deriv(Γ::AbstractMetricGraph, σ::QuantumGraphEigen, coef::Vector)
    Λ_sqrt = sqrt.(σ.Λ);
    u_Q_deriv = Vector{Function}(undef,ne(Γ.G));
    for j = 1:ne(Γ.G)
        function u_Q_deriv_e(x)
            Ae_Cos_x=[-σ.A[j,q]*sin(Λ_sqrt[q]*x)*Λ_sqrt[q] for q=1:σ.Q];
            Be_Sin_x=[σ.B[j,q]*cos(Λ_sqrt[q]*x)*Λ_sqrt[q] for q=1:σ.Q];
            return coef'*(Ae_Cos_x+Be_Sin_x)
        end
        u_Q_deriv[j] = u_Q_deriv_e;
    end
    return u_Q
end

"""
    spectral_solver(TP::EllipticTestProblem, K::Int; filon=false)

Compute coefficents of the spectral Galerkin solution of 'TP' with eigenfunction basis of size 'K'.

If filon=true, the more economic filon-quadrature is applied instead of QuadGK.
"""
function spectral_solver(TP::EllipticTestProblem, K::Int; filon=false)
    σ = eigen_quantum(TP.Γ, K=K);
    if filon == false
        f = projection_coefficients(TP.Γ, σ, TP.rhs)
    else
        f = projection_coefficients_filon(TP.Γ, σ, TP.rhs, 15*K)
    end
    return Diagonal((σ.Λ + TP.pot*ones(σ.Q)).^(-1)) * f
end

"""
    spectral_solver(TP::TPGeneralizedHeat, T::Number, K::Int; filon=false)

Compute coefficents of the spectral Galerkin solution of 'TP' at time 'T' with eigenfunction basis of size 'K'.

If filon=true, the more economic filon-quadrature is applied instead of QuadGK.
"""
function spectral_solver(TP::TPGeneralizedHeat, T::Number, K::Int; filon=false)
    σ = eigen_quantum(TP.Γ, K=K);
    if filon == false
        coef = projection_coefficients(TP.Γ, σ, TP.u0)
        return diagm(exp.(-T*σ.Λ))*coef
    else
        coef = projection_coefficients_filon(TP.Γ, σ, TP.u0, 15*K)
        return diagm(exp.(-T*σ.Λ))*coef
    end
end


"""
    trunc_error(TP::EllipticTestProblem, σ::QuantumGraphEigen, coef::Vector)

"""
function trunc_error(TP::EllipticTestProblem, σ::QuantumGraphEigen, coef::Vector)
    spec_sol = spectral_solution(TP.Γ, σ, coef);
    if typeof(TP.Γ) == EquilateralMetricGraph
        ℓ_vec = TP.Γ.ℓ * ones(Int, ne(TP.Γ.G))
    else 
        ℓ_vec = TP.Γ.ℓ_vec
    end
    trunc_err = 0;
    for j = 1:ne(TP.Γ.G)
        trunc_diff(x) = (TP.u[j](x)-spec_sol[j](x))^2;
        trunc_err += quadgk(trunc_diff,0,ℓ_vec[j])[1]
    end
    return sqrt(trunc_err)
end


"""
    L2_norm_spectral(coef::Vector)

Compute L^2 Norm of spectral solution ``u_Q = \\sum_{q \\le Q} 'coef'_q \\phi_q``

"""
function L2_norm_spectral(coef::Vector)
    return sqrt(coef'*coef)
end

"""
    H1_seminorm_spectral(σ::QuantumGraphEigen, coef::Vector)

Computes ``H^1(\\Gamma)`` semimorm of spectral solution ``u_Q = \\sum_{q \\le Q} 'coef'_q \\phi_q`` with eigenvalues 'σ.Λ'.

"""
function H1_seminorm_spectral(σ::QuantumGraphEigen, coef::Vector)
    return sqrt(coef'*Diagonal(σ.Λ)*coef)
end

"""
    H1_seminorm_spectral(Λ::Vector, coef::Vector)

Computes `H^1(\\Gamma)`` of spectral solution ``u_Q = \\sum_{q \\le Q} 'coef'_q \\phi_q`` with eigenvalues 'Λ'.

"""
function H1_seminorm_spectral(Λ::Vector, coef::Vector)
    return sqrt(coef'*Diagonal(Λ)*coef)
end