# Iterative Solution Method in PyFEMP
In the kernel of PyFEMP the solution of a nonlinear system of equations is computed. This solution is obtained by a Newton-Raphson iterative algorithm.

Generally it is assumed that any system of equations to solve, is given in a residual form defined as a functions of a vector of unknowns

```math
\boldsymbol{R}(\boldsymbol{d}) = \boldsymbol{0}.
```

This system of equations by default is solved iteratively, based on the first order Taylor expansion 

```math
\boldsymbol{R}(\boldsymbol{d}) + \dfrac{\text{d} \boldsymbol{R}}{\text{d} \boldsymbol{d}} \, \Delta \boldsymbol{d} = \boldsymbol{0},
```
introducing an incremental change $`\Delta \boldsymbol{d}`$.

In terms of the finite elment discretization, this residual is composed of element contributions (i.e. the assembly of elememt vectors) assembelec to the global vector $`\boldsymbol{P}`$ and a global load vector build from conventional nodal loads (natural boundray conditions/ Neumann boundary conditions) $`\boldsymbol{F}`$.

```math
\boldsymbol{R}(\boldsymbol{d}) = \boldsymbol{F} - \boldsymbol{P}(\boldsymbol{d}) = \boldsymbol{0}.
```

Applying this definition of the residual function to the iterative scheme

```math
\boldsymbol{F} - \boldsymbol{P}(\boldsymbol{d}) +
\dfrac{\text{d} \boldsymbol{R}}{\text{d} \boldsymbol{d}} \, \Delta \boldsymbol{d}  = \boldsymbol{0},
```

we obtain the actually solved linear system

```math
\boldsymbol{K} \, \Delta \boldsymbol{d}  = \boldsymbol{F} - \boldsymbol{P}(\boldsymbol{d}),
```

introducing the consistent tangent

```math
\boldsymbol{K} = - \dfrac{\text{d} \boldsymbol{R}}{\text{d} \boldsymbol{d}}.
```



