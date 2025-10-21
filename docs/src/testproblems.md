# Test Problems

Different type of test problems of abstract type TestProblem can be constructed. 

The following types are available:

- [`EllipticTestProblem`](@ref)
- [`TPGeneralizedHeat`](@ref)
- [`TPReactionDiffusion`](@ref)

Depending on the specific problem, a right-hand side, reaction-term and/or initial conditions must be set. 
If available, the exact solution and its derivative can be indicated. Otherwise, 'Nothing' must be specified. 

## Elliptic Test Problems

The elliptic equation 

```math
\mathcal{H}( u(x)) + \nu u(x) = f(x)
```

is represented by the TestProblem type 

```@docs
EllipticTestProblem
```

The following exemplary test problems are predefined:

```@docs
TestProblem242
TestProblem243
TestProblem711
```

## Parabolic Test Problems

### Generalized Heat Equation

The generalized heat equation

```math
\frac{\partial u}{\partial t}(x,t) + \mathcal{H}( u(x,t) = f(x,t))
```

on ``\Gamma`` subject to Neumann-Kirchhoff conditions and initial condition `` u^0 `` is represented by the TestProblem type 

```@docs
TPGeneralizedHeat
```

The following exemplary test problems are predefined:

```@docs
TestProblem244
``` 

```@docs
TestProblem245
``` 

```@docs
TestProblem721
``` 

### Reaction-Diffusion Equation

```@docs
TPReactionDiffusion
```
