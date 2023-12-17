# ... and its spectrum

All eigenvalues ``\lambda < \frac{K \pi}{\ell}^2`` of the equilateral graph ``\Gamma`` can be computed with [`eigvals_quantum`](@ref). By default, K=3 is applied. 

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
using Plots # hide
Γ = metric_star_graph() # hide
eigvals_quantum(Γ) 
```

An eigenfunction basis with all eigenfunctions ``\phi_{\lambda}``, ``\lambda < \frac{K \pi}{\ell}^2`` can be constructed via [`eigen_quantum`](@ref)

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
using Plots # hide
Γ = metric_star_graph() # hide
σ = eigen_quantum(Γ) 
```

Eigenvalues and eigenfunctions are always returned in ascending order. 
The function  allows to explicitly construct a specific eigenfunction:

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
using Plots # hide
Γ = metric_star_graph() # hide
σ = eigen_quantum(Γ) # hide
ϕ_q = eigenfunction(Γ, σ, 5)
```

It can be vizualized using [`plot_function_3d`](@ref)

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
using Plots # hide
Γ = metric_star_graph() # hide
σ = eigen_quantum(Γ) # hide
ϕ_q = eigenfunction(Γ, σ, 5) # hide
plot_function_3d(Γ, ϕ_q)
```