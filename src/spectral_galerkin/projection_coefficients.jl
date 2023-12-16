# © Anna Weller, University of Cologne, 2023

using QuadGK
using Graphs

### Basic implementation with QuadGK

"""
    projection_coefficients(Γ::MetricGraph, σ::QuantumGraphEigen, u::Vector{Function})

Compute projection coefficients 'coefs' of orthogonal projection of 'u' with standard QuadGK quadrature.
"""
function projection_coefficients(Γ::MetricGraph, σ::QuantumGraphEigen, u::Vector{Function})
    coefs = zeros(σ.Q)
    # compute projection coefficients for all lambda
    for q = 1:σ.Q
        λ = σ.Λ[q]; √λ = sqrt(λ);
        coef_lambda=0
        for j = 1:ne(Γ.G) 
            A_e = σ.A[j,q]; B_e = σ.B[j,q]
            ϕ_e(x) = A_e*cos(√λ*x)+B_e*sin(√λ*x)
            # coef_func(x)=f_to_approx[edge](x)*eigfunc(x)
            coef_func(x) = u[j](x)*ϕ_e(x)
            integral,_ = quadgk(coef_func, 0, Γ.ℓ_vec[j])
            coef_lambda += integral
        end
        coefs[q] = coef_lambda
    end
    return coefs
end

"""
    projection_coefficients(Γ::EquilateralMetricGraph, σ::QuantumGraphEigen, u::Vector{Function})

Equilateral version.
"""
function projection_coefficients(Γ::EquilateralMetricGraph, σ::QuantumGraphEigen, u::Vector{Function})
    coefs = zeros(σ.Q)
    # compute projection coefficients for all lambda
    for q = 1:σ.Q
        λ = σ.Λ[q]; √λ = sqrt(λ);
        coef_lambda=0
        for j = 1:ne(Γ.G) 
            A_e = σ.A[j,q]; B_e = σ.B[j,q]
            ϕ_e(x) = A_e*cos(√λ*x)+B_e*sin(√λ*x)
            # coef_func(x)=f_to_approx[edge](x)*eigfunc(x)
            coef_func(x) = u[j](x)*ϕ_e(x)
            integral,_ = quadgk(coef_func, 0, Γ.ℓ)
            coef_lambda += integral
        end
        coefs[q] = coef_lambda
    end
    return coefs
end


## Matrix-free computation via Filon-Quadtrature 


"""
    projection_coefficients_filon(Γ::EquilateralMetricGraph, σ::QuantumGraphEigen, u::Vector{Function}, quad_nodes::Int)

Compute projection coefficients of 'u' on 'Γ' for eigenfunction in 'σ' with 'quad_nodes' quadrature nodes.

Uses (matrix-free) Filon-Quadrature as described in \\[W\\], Section 4.3.2.

""" 
function projection_coefficients_filon(Γ::EquilateralMetricGraph, σ::QuantumGraphEigen, u::Vector{Function}, quad_nodes::Int)
    h_kap = Γ.ℓ/quad_nodes; N_kap = Int(round(Γ.ℓ/h_kap))
    coefs = zeros(σ.Q)
    fill_proj_coefs!(coefs,σ,u,h_kap,N_kap)
    return coefs
end

"""
    fill_proj_coefs!(coefs::Vector, σ::QuantumGraphEigen, u::Vector{Function}, h_kap::Number, N_kap::Int)

"""
function fill_proj_coefs!(coefs::Vector, σ::QuantumGraphEigen, u::Vector{Function}, h_kap::Number, N_kap::Int)
    for j=1:length(u)
        vec_x2,vec_x,vec_1 = function_coef_evaluation(h_kap, N_kap, u, j);
        if iszero(vec_x2) && iszero(vec_x) && iszero(vec_1)
            continue
        else
        #
            coefs[1]+=σ.A[j,1]*mom_integralsf_evals_q0(h_kap, N_kap, vec_x2, vec_x, vec_1, 0, 0)
            for q = 2:σ.Q
                λ = σ.Λ[q]; λ_sqrt = sqrt(λ)
                if isapprox(σ.A[j,q],0,atol=10^(-14)) && isapprox(σ.B[j,q],0,atol=10^(-14))
                    continue
                else
                    coefs[q]+=integrals_q_j(h_kap, N_kap, vec_x2, vec_x, vec_1, λ, λ_sqrt, σ.A[j,q], σ.B[j,q], q)
                end    
            end
        end
    end
end

"""
    fill_vecs!(vec_x2::Vector, vec_x::Vector, vec_1::Vector, h_kap::Number, N_kap::Int, boldsymbolu_e, coefs_alpha, count)

"""
function fill_vecs!(vec_x2::Vector, vec_x::Vector, vec_1::Vector, h_kap::Number, N_kap::Int, boldsymbolu_e, coefs_alpha, count)
    for kap = 1:2:N_kap-1
        f_evals_kap = boldsymbolu_e[kap:kap+2];
        vec_x2[count] = coefs_alpha'*f_evals_kap
        vec_x[count] = [-(kap+0.5)/(h_kap), 2*kap/h_kap, -(kap-0.5)/(h_kap)]'*f_evals_kap
        vec_1[count] = ([(kap^2+kap)/2,-(kap-1)*(kap+1),(kap^2-kap)/2])'*f_evals_kap
        count += 1
    end
end

"""
    function_coef_evaluation(h_kap::Number, N_kap::Int, u::Vector{Function}, j::Int)

Compute function evaluations of 'u'_'j' on quadrature nodes x_kap = kap*'h_kap', kap=1,…'N_kap'
"""
function function_coef_evaluation(h_kap::Number, N_kap::Int, u::Vector{Function}, j::Int)
    boldsymbolu_e = u[j].([count*h_kap for count = 0:N_kap]);
    if iszero(boldsymbolu_e)
        vec_x2 = zeros(Int(N_kap/2)); vec_x = zeros(Int(N_kap/2)); vec_1 = zeros(Int(N_kap/2))
        return vec_x2, vec_x, vec_1
    else
        vec_x2 = zeros(Int(N_kap/2)); vec_x = zeros(Int(N_kap/2)); vec_1 = zeros(Int(N_kap/2))
        fill_vecs!(vec_x2, vec_x, vec_1, h_kap, N_kap, boldsymbolu_e, [1/(2*h_kap^2),-1/(h_kap^2),1/(2*h_kap^2)], 1)
        return vec_x2, vec_x, vec_1
    end
end

"""
    mom_integralsf_evals_q0(h_kap::Number, N_kap::Int, vec_x2::Vector, vec_x::Vector, vec_1::Vector)

Compute moment integrals.
"""
function mom_integralsf_evals_q0(h_kap::Number, N_kap::Int, vec_x2::Vector, vec_x::Vector, vec_1::Vector, λ::Number, λ_sqrt::Number)
    summe = 0; count = 1
    for kap = 1:2:N_kap-1
        # I_1,iota
        I_12_t = 1/3*((h_kap*(kap+1))^3-(h_kap*(kap-1))^3)
        I_11_t = 1/2*((h_kap*(kap+1))^2-(h_kap*(kap-1))^2)  
        I_10_t = 2*h_kap
        summe+=I_12_t*vec_x2[count]+I_11_t*vec_x[count]+I_10_t*vec_1[count]
        count+=1
    end
    return summe
end

"""
    integrals_q_j(h_kap::Number, N_kap::Int, vec_x2::Vector, vec_x::Vector, vec_1::Vector, λ, λ_sqrt, Ajq, Bjq, q)

"""
function integrals_q_j(h_kap::Number, N_kap::Int, vec_x2::Vector, vec_x::Vector, vec_1::Vector, λ, λ_sqrt, Ajq, Bjq, q)
    if isapprox(Ajq,0,atol=10^(-14))
        sum2 = 0;
        # compute only sum2 
        count = 1;
        for kap = 1:2:N_kap-1
            c0 = (1/(λ^(3/2))); c1 = (h_kap^2*(kap-1)^2*λ-2); c2 = (h_kap^2*(kap+1)^2*λ-2); c3 = (2/λ); c4 = (2/λ_sqrt)
            arg1 = h_kap*(kap-1)*λ_sqrt; arg2 = h_kap*(kap+1)*λ_sqrt; arg3 = h_kap*kap*λ_sqrt; arg4 = h_kap*λ_sqrt
            cosarg1 = cos(arg1); cosarg2 = cos(arg2); cosarg3 = cos(arg3); cosarg4 = cos(arg4)
            sinarg1 = sin(arg1); sinarg2 = sin(arg2); sinarg3 = sin(arg3); sinarg4 = sin(arg4)
            I_22_t = c0*(c1*cosarg1-c2*cosarg2-2*arg1*sinarg1+2*arg2*sinarg2)
            I_21_t = c3*(sinarg4*(arg3*sinarg3+cosarg3)-arg4*cosarg4*cosarg3)
            I_20_t = c4*(sinarg4*sinarg3)
            sum2 += I_22_t*vec_x2[count]+I_21_t*vec_x[count]+I_20_t*vec_1[count]
            count += 1;
        end
        #proj_coefs_edge[q]=Bjq*sum2
        return Bjq*sum2
    elseif isapprox(Bjq,0,atol=10^(-14))
        sum1 = 0;
        count = 1;
        for kap = 1:2:N_kap-1
            c0 = (1/(λ^(3/2))); c1 = (h_kap^2*(kap-1)^2*λ-2); c2 = (h_kap^2*(kap+1)^2*λ-2); c3 = (2/λ); c4 = (2/λ_sqrt)
            arg1 = h_kap*(kap-1)*λ_sqrt; arg2 = h_kap*(kap+1)*λ_sqrt; arg3 = h_kap*kap*λ_sqrt; arg4 = h_kap*λ_sqrt
            cosarg1 = cos(arg1); cosarg2 = cos(arg2); cosarg3 = cos(arg3); cosarg4 = cos(arg4)
            sinarg1 = sin(arg1); sinarg2 = sin(arg2); sinarg3 = sin(arg3); sinarg4 = sin(arg4)
    
            I_12_t = c0*(-c1*sinarg1+ c2*sinarg2- 2*arg1*cosarg1+ 2*arg2*cosarg2)
            I_11_t = c3*(arg3*sinarg4*cosarg3+sinarg3*(arg4*cosarg4-sinarg4))     
            I_10_t = c4*(sinarg4*cosarg3)
            sum1 += I_12_t*vec_x2[count]+I_11_t*vec_x[count]+I_10_t*vec_1[count]
            count += 1;
        end
        #proj_coefs_edge[q]=Ajq*sum1
        return Ajq*sum1
    else
        sum1 = 0; sum2 = 0
        count = 1;
        for kap=1:2:N_kap-1
            c0 = (1/(λ^(3/2))); c1 = (h_kap^2*(kap-1)^2*λ-2); c2 = (h_kap^2*(kap+1)^2*λ-2); c3 = (2/λ); c4 = (2/λ_sqrt)
            arg1 = h_kap*(kap-1)*λ_sqrt; arg2 = h_kap*(kap+1)*λ_sqrt; arg3 = h_kap*kap*λ_sqrt; arg4 = h_kap*λ_sqrt
            cosarg1 = cos(arg1); cosarg2 = cos(arg2); cosarg3 = cos(arg3); cosarg4 = cos(arg4)
            sinarg1 = sin(arg1); sinarg2 = sin(arg2); sinarg3 = sin(arg3); sinarg4 = sin(arg4)
    
            I_12_t = c0*(-c1*sinarg1+ c2*sinarg2- 2*arg1*cosarg1+ 2*arg2*cosarg2)
            I_11_t = c3*(arg3*sinarg4*cosarg3+sinarg3*(arg4*cosarg4-sinarg4))     
            I_10_t = c4*(sinarg4*cosarg3)
            sum1 += I_12_t*vec_x2[count]+I_11_t*vec_x[count]+I_10_t*vec_1[count]
    
            I_22_t = c0*(c1*cosarg1 - c2*cosarg2 - 2*arg1*sinarg1 + 2*arg2*sinarg2)
            I_21_t = c3*(sinarg4*(arg3*sinarg3+cosarg3)-arg4*cosarg4*cosarg3)
            I_20_t = c4*(sinarg4*sinarg3)
            sum2 += I_22_t*vec_x2[count]+I_21_t*vec_x[count]+I_20_t*vec_1[count]
            count += 1
        end
        #proj_coefs_edge[q]=Ajq*sum1+Bjq*sum2
        return Ajq*sum1+Bjq*sum2
    end
end
