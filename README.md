# MeGraPDE

MeGraPDE stands for MetricGraphPDEs and implements numerical methods for the solution of partial differential equations (PDEs) on metric graphs. 

I have developed the MeGraPDE.jl package in connection with my Ph.D thesis cite at the University of Cologne \[W\]. 

Installation:

Please see the Documentation for usage instructions.

Among others, the package includes

- construction of a variety of exemplary metric graphs and test problems
- discretization via extended graphs
- finite element solver for PDEs on metric graphs, e.g. in combination with a multigrid approach
- computation of quantum graph eigenvalues and eigenfunctions
- spectral Galerkin solver for PDEs on metric graphs, e.g. in combination with a filon-quadrature

The finite element discretization via extended graphs is implemented based on the original work \[AB\]. 
The computation of equilateral quantum graph eigenvalues is based on an idea originally proposed by von Below \[B\].
The remaining methods and the related theory have been derived for \[W\] and are discussed therein.

The package is under continuous development. 

\[W\] Anna Weller, Numerical Methods for Parabolic Partial Differential Equations on Metric Graphs, PhD thesis at the University of Cologne, in preparation.

\[AB\] Mario Arioli, Michele Benzi, A finite element method for quantum graphs, IMA Journal of Numerical Analysis, Volume 38, Issue 3, July 2018, Pages 1119–1163.

\[B\] Joachim von Below, A characteristic equation associated to an eigenvalue problem on c2-networks. Linear Algebra and its Applications, 71:309–325, 1985.

Copyright (c) 2023 Anna Weller (University of Cologne)

