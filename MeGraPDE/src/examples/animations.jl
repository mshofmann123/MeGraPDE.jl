# © Anna Weller, University of Cologne, 2023


using Plots, Measures
using Graphs, LinearAlgebra

function animate_diffusion(Γ::EquilateralMetricGraph, u0::Vector{Function}, T::Number; α = 1)
    dt = T /50
    σ = eigen_quantum(Γ, K=20)
    coefs = projection_coefficients_filon(Γ, σ, u0, 500)
    anim = @animate for t in 0:dt:T
        coef_t = exp(-t*Diagonal(σ.Λ.^(α)))*coefs
        u_Q=spectral_solution(Γ, σ, coef_t)
        plot_function_3d(Γ, u_Q, size=(200,200), lw=3, grid_off=true)
    end
    gif(anim, fps=5)
end

