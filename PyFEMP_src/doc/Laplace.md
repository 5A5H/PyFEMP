# Steady State Heat Conduction - The Laplace Problem

>> Please recognize that this is intended to be only a short summary to illustrate the implementation of elements in PyFEMP. It might be slopy written and incomplete in terms of theoretical aspects.

The steady state heat conductivity problem is given in the strong form as
```math
\text{div}(\boldsymbol{q}) = \boldsymbol{0},
```
with the heat flux $`\boldsymbol{q} = - \alpha_{q} \,\text{grad} \theta`$. Hereby $`\alpha_{q} \geq 0`$ denotes the heat conductivity parameter and $`\theta`$ the absolute temperature.

To obtain an approximate solution using the finite element method, we build a weak form using the arbitrary test function $`\delta\theta`$, arriving at

```math
G = \int_{\mathcal B} \, \boldsymbol{q} \cdot \text{grad}\, \delta \theta \, \text{d}v - \int_{\partial {\mathcal B}^{q}} \bar{\boldsymbol{q}} \cdot \boldsymbold{n} \, \delta \theta \,\text{d}a
```