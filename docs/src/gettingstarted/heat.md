# Heat Equation

Finally, the initial boundary value problem for the heat equation

```math
\frac{\partial}{\partial t} u(x,t) + \frac{\partial^2}{\partial u^2} u(x,t) = 0 \quad (\ast)
```

on a metric graph `Γ` under Neumann-Kirchhoff conditions is approximated.

Consider a lollipop graph that can be constructed using the predefined method [`metric_lollipop_graph`](@ref)

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
Γ = metric_lollipop_graph()
plot_graph_3d(Γ)
```

As initial condition, we choose a model initial condition that has compact support on randomly chosen, edge of `Γ ` and is zero elsewhere. A routine to assemble this initial condition is implemented in 

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
Γ = metric_lollipop_graph() # hide
u0 = model_initial_condition(Γ) 
```

The solution of ``(\ast)`` can be simulated for ``t \in [0,T]`` by calling

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
Γ = metric_lollipop_graph() # hide
u0 = model_initial_condition(Γ) # hide
T = 1
animate_diffusion(Γ, u0, T)
```

Lets go for fractional diffusion on a tree!

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
Γ = metric_tree_graph()
u0 = model_initial_condition(Γ)
T = 3
animate_diffusion(Γ, u0, T, α = 0.1)
```



 