# Base

The abstract type AbstractMetricGraph comprises the following two types representing non-equilateral and equilateral metric graphs.

```@docs
    MetricGraph
    EquilateralMetricGraph
```

In the field 'coords', coordinates can be specified that allow plotting the graph in a 3d grid.

The function metric_graph assembles a metric graph by specifying a combinatorial graph `G` and edge lengths. Coordinates can be set optionally, but are false per default.

```@docs
    metric_graph
```

It is possible to extended every possible graph to a metric graph by assigning edge length. 

Consider for example the Barabasi-Albert graph constructed with [Graphs.jl](https://github.com/JuliaGraphs/Graphs.

```@example
using Graphs
G = barabasi_albert(100,3)
```

You can now either assign an equilateral edge length represented by one number or a vector l_vec with edge lengths to create a metric graph

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
using Graphs #hide
G = barabasi_albert(100,3) #hide
ℓ = 1
Γ = metric_graph(G, ℓ)
```

```@example
push!(LOAD_PATH,"../src/") # hide
using ..MeGraPDE # hide
using Graphs #hide
G = barabasi_albert(100,3) #hide
ℓ_vec = rand(1:5,ne(G))
Γ = metric_graph(G, ℓ_vec)
```

Some minor functions are implemented to quickly access properties of Γ. This list is by far not complete and will be expanded by other frequently used functions. 

```@docs
edge_length
vol
```



