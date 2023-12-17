# MeGraPDE

MeGraPDE stands for MetricGraphPDEs and implements numerical methods for the solution of partial differential equations (PDEs) on metric graphs. I have developed the MeGraPDE.jl package in connection with my Ph.D thesis at the University of Cologne \[W\]. 

## Documentation

Full documentation and first usage examples are available at https://annaweller.github.io/MeGraPDE.jl/

## Info

Among others, the package includes

- construction of a variety of exemplary metric graphs and test problems
- discretization via extended graphs
- finite element solver for PDEs on metric graphs, e.g. in combination with a multigrid approach
- computation of quantum graph eigenvalues and eigenfunctions
- spectral Galerkin solver for PDEs on metric graphs, e.g. in combination with a filon-quadrature

The package relies on the methods from [Graphs.jl](https://github.com/JuliaGraphs/Graphs.jl) for combinatorial graphs.

The finite element discretization via extended graphs is implemented based on the original work \[AB\]. 
The computation of equilateral quantum graph eigenvalues is based on an idea originally proposed by von Below \[B\].
The remaining methods and the related theory have been derived for \[W\] and are discussed therein.

The package is under continuous development. 

\[W\] Anna Weller, Numerical Methods for Parabolic Partial Differential Equations on Metric Graphs, PhD thesis at the University of Cologne, in preparation.

\[AB\] Mario Arioli, Michele Benzi, A finite element method for quantum graphs, IMA Journal of Numerical Analysis, Volume 38, Issue 3, July 2018, Pages 1119–1163.

\[B\] Joachim von Below, A characteristic equation associated to an eigenvalue problem on c2-networks. Linear Algebra and its Applications, 71:309–325, 1985.

Copyright (c) 2023 Anna Weller (University of Cologne)

## Installation

The package can be added by specifying the URL to the Git repository. In your `julia` terminal, enter the following commands


```julia
julia> using Pkg
julia> Pkg.add(url="https://github.com/AnnaWeller/MeGraPDE.jl");
```

You are all set. The package can now be activated with the command 

```@repl
using MeGraPDE
```

