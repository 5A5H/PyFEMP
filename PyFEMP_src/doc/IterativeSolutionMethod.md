# Iterative Solution Method in PyFEMP
In the kernel of PyFEMP the solution of a nonlinear system of equations is computed. This solution is obtained by a Newton-Raphson iterative algorithm.

Generally it is assumed that any system of equations to solve, is given in a residual form defined as a functions of a vector of unknowns

```math
\boldsymbol{R}(\boldsymbol{d}) = \boldsymbol{0}.
```

This system of equations by default is solved iteratively, based on the first order taylor expansion 

```math
\boldsymbol{R}(\boldsymbol{d}) + \Delta \boldsymbol{d} \, \dfrac{\text{d} \boldsymbol{R}}{\text{d} \boldsymbol{d}} = \boldsymbol{0},
```
introducing an incremental change `$\boldsymbol{d}$`.