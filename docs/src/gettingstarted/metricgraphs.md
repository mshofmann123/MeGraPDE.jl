# Your first metric graph

Before you start, go to the installation section and activate `MeGraPDE` in the current session with the command

```@example
using MeGraPDE
```

## Creating a metric graph

Let us first create a combinatorial star graph `G` with 5 vertices and 4 edges using [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl)

```@example
using Graphs
G = star_graph(5)
```

In order to extend `G` to an equilateral metric graph, we define the edge length `ℓ`

```@example
ℓ = pi + pi/2
```

`G` can now be represented as metric graph `Γ` by applying the function [`metric_graph`](@ref)

```@example
push!(LOAD_PATH, "../src/") # hide
using ..MeGraPDE # hide
using Graphs # hide
G = star_graph(5) # hide
ℓ = pi + pi/2 #hide
Γ = metric_graph(G, ℓ)
```

For a small example like the star graph, vertex coordinates can be assigned that will later allow to visualize `Γ` in 3d.

```@example
push!(LOAD_PATH, "../src/") # hide
using ..MeGraPDE # hide
using Graphs # hide
G = star_graph(5) # hide
ℓ = pi + pi/2 #hide
coords = [[0, 0], [ℓ, 0], [-ℓ, 0], [0, ℓ], [0, -ℓ]]
```

The function [`metric_graph`](@ref) takes the optional input `vertex_coords` to specify the vertex coordinates.

```@example
push!(LOAD_PATH, "../src/") # hide
using ..MeGraPDE # hide
using Graphs # hide
G = star_graph(5) # hide
ℓ = pi + pi/2 #hide
coords = [
    [0, 0], #hide
    [ℓ, 0], #hide
    [-ℓ, 0], #hide
    [0, ℓ], #hide
    [0, -ℓ],
] # hide
Γ = metric_graph(G, ℓ; vertex_coords=coords)
```

We may now plot `Γ` using [`plot_graph_3d`](@ref)

```@example
push!(LOAD_PATH, "../src/") # hide
using ..MeGraPDE # hide
using Plots # hide
Γ = metric_star_graph() # hide
plot_graph_3d(Γ)
```

!!! note
    

The previous example graph can be assembled using the constructor [`metric_star_graph`](@ref)
and indicating the desired edge length as `metric_star_graph(ℓ = pi + pi/2)`.
Several other example graphs are implemented.

## Functions on metric graphs

A function `u` on a metric graph is represented by a vector of functions `u_e`, specifying `u` on each edge `e`.

```@example
u = [x -> -3*sin(x), x -> sin(x), x -> sin(x), x -> sin(x)]
```

If vertex coordinates are assigned to `Γ`, a function can be plotted on `Γ` with

```@example
push!(LOAD_PATH, "../src/") # hide
using ..MeGraPDE # hide
using Graphs # hide
G = star_graph(5) # hide
ℓ = pi + pi/2 #hide
coords = [
    [0, 0], #hide
    [ℓ, 0], #hide
    [-ℓ, 0], #hide
    [0, ℓ], #hide
    [0, -ℓ],
] # hide
Γ = metric_graph(G, ℓ; vertex_coords=coords) # hide
u = [
    x -> -3*sin(x), # hide
    x -> sin(x),  # hide
    x -> sin(x), # hide
    x -> sin(x), # hide
] # hide
plot_function_3d(Γ, u)
```
