# Steady State Heat Conduction - The Laplace Problem

>> Please recognize that this is intended to be only a short summary to illustrate the implementation of elements in PyFEMP. It might be slopy written and incomplete in terms of theoretical aspects.

The steady state heat conductivity problem is given in the strong form as
```math
\text{div}(\boldsymbol{q}) = \boldsymbol{0},
```
with the heat flux $`\boldsymbol{q} = - \alpha_{q} \,\text{grad} \theta`$. Hereby $`\alpha_{q} \geq 0`$ denotes the heat conductivity parameter and $`\theta`$ the absolute temperature.

To obtain an approximate solution using the finite element method, we build a weak form using the arbitrary test function $`\delta\theta`$, arriving at

```math
G = 
\int_{\mathcal B} \, \boldsymbol{q} \cdot \text{grad} \, \delta \theta \, \text{d}v - \int_{\partial {\mathcal B}^{q}} \, \bar{\boldsymbol{q}} \cdot \boldsymbol{n} \, \delta \theta \,\text{d}a
```
where $`\bar{\boldsymbol{q}} \cdot \boldsymbol{n}`$ is the applied, heat flux on the free boundary (no essential/dirichlet condition) $`\partial {\mathcal B}^{q}`$. Remark: A positive value ($`\bar{\boldsymbol{q}} \cdot \boldsymbol{n} > 0`$) does indicate a heatflux out of the body, due to the definition of the outward pointing normal $`\boldsymbol{n}`$.

# PyFEMP - Finite elements
A finite element for usage in PyFEMP is a single python file which defines three functions:
 * `Elmt_Init()` which tells PyFEMP about the Layout of the element (e.g. names/numbers of nodal degrees of freedom, shape, material parameters)
 * `Elmt_KS(XL, UL, Hn, Ht, Mat, dt)` for computing element vector and matrix
 * `Elmt_Post(XL, UL, Hn, Ht, Mat, dt, PostName)` for computing postprocessing vectors containing a value for each node

The element file must `import numpy as np` in the beginning.

There is no further requirements or restrictions to a element file. You may define addditional function, put comments or else. Especially feel free to add some `print()` statements if you are curious about something.

# A finite element for the steady state Laplace problem
In the following we will discuss the essential functions to be implemented for computing the laplace problem described above.

# The `Elmt_Init()`
This function is called by PyFEMP in the beginning and tells a simulation about the element properties.
```python
def Elmt_Init():
    NoElementDim = 2                    # Number of dimensions
    NoElementNodes = 4                  # Number of nodes on the element (here: 4-Node Quadirlateral)
    NoElementHistory = 0                # Number of history variables
    ElementDofNames = ["T"]             # Element degrees of freedom specified by a string for each scalar
    ElementMaterialNames = ["alpha_q"]  # Element material parameters specified by a string for each scalar
    ElementPostNames = ["T"]            # Element postprocessing parameters
    return NoElementDim, NoElementNodes, ElementDofNames, NoElementHistory, ElementMaterialNames, ElementPostNames
```

In this example we will build a 4-node quadiralteral element (`NoElementNodes = 4`) for a 2D (`NoElementDim = 2`) analysis. As we have a steady state problem, we don't have any requirement for element history variables (`NoElementHistory = 0`). 
In PyFEMP it is assumed that each node of an element carries the same number of scalar degrees of freedom. In the `Elmt_Init()` function these are specified by a string, which will also be used during the analysis (see: XX).
Here each element node has only one scalar degree of freedom, which is the absolute temperature (`ElementDofNames = ["T"]`). 
Also the list of material parameters is rather short on this example, where we only need the heat conductivity parameter (`ElementMaterialNames = ["alpha_q"]`).
Postprocessing fields are introduced by a string name per scalar as well, but we'll not go into details here.

# The `Elmt_KS(XL, UL, Hn, Ht, Mat, dt)`

## Interface of the `Elmt_KS` routine

This is the main routine of the finite element, called during the solution procedure. The input is standartized as

 * `XL = np.array([X1x, X1y, X2x, X2y, X3x, X3y, X4x, X4y])` list of nodal coordinates. Here we have 4 nodes with an x, and y coordinate for each node.
 Hence its size is `NoElementDim * NoElementNodes`.
 
 * `UL = np.array([T1, T2, T3, T4])` list of current dof value. That means, the element gets the current temperature within the iterative scheme as $`\theta = \theta^{t} + \Delta \theta`$ where $`\theta `$ is the temperature in the element, $`\theta^{t}`$ is the value of the temperature at the beginning of the time step (where we don't have equilibrium yet) and $`\Delta \theta`$ is the increment towards equilibrium in the current time step. Its size is `len(ElementDofNames) * NoElementNodes`.
 
 * `Hn`, `Ht` are vectors with the element history variables. `Hn` contain the values at the beginning of the time/load step and `Ht` is its current counter part. Both are of equal size, `NoElementHistory`.
 
 * `Mat` is a vector with the material parameters, in our case with one entry for the conductivity. Generally its size is `len(ElementMaterialNames)`.
 
 * `dt` is a scalar, being the time increment from the beginning of the time step to its end.


 >> Note that the lists `XL` and `UL` are always ordered nodewise. For a coupled element with a displacement vector and a temperature as nodal degrees of freedom we would have: `UL = np.array([U1x, U1y, T1, ... , Unx, Uny, Tn])` for `n`- nodes. 

Using the standard element procedure the element vector and matrix are computed, and returned. Recognize that the element routine does allocate its vector and matrix as numpy arrays, within the `Elmt_KS` routine, i.e.
```python
    # initialize element vector /matrix
    r_e = np.zeros(NoElmtNodes*NoNodalDOF)
    k_e = np.zeros((NoElmtNodes*NoNodalDOF, NoElmtNodes*NoNodalDOF))
```

Hereby the PyFEMP notation demands that the weak form to be solved
is directly implemented to the element vector (without changing the sign).
That means specifically, the element vector $`\boldsymbol{r}_{e}`$ is defined as:

```math
\boldsymbol{r}_{e} = \dfrac{\text{d}}{\text{d} \delta \boldsymbol{d}_{e}} \, G
```
with the vector of nodal tests $`\delta \boldsymbol{d}_{e}`$.

Furthermore the element matrix is supposed to be the total derivative of the element vector, w.r.t. the vector of nodal dofs $`\boldsymbol{d}_{e}`$. In the routine this is represented by the `UL` input vector.

```math
\boldsymbol{k}_{e} = \dfrac{\text{d}}{\text{d} \boldsymbol{d}_{e}} \, \boldsymbol{r}_{e}
```


## Complete code of the `Elmt_KS` routine
```python
def Elmt_KS(XL, UL, Hn, Ht, Mat, dt):
    '''
    This function returns element vector and matrix.
    Generally the element matrix is the straight derivative
    of the element vector with respect to the nodal degrees of freedoms,
    in the order of the UL field.
    Both need to be returned as numpy arrays, with dimensions: 
    element vector r_e : number of nodal degrees of freedom times number of nodes
    element matrix k_e : NoNodes * NoNodalDOF ,NoNodes * NoNodalDOF
    '''
    verbose = False  # True;
    if verbose: print('XI :',XL)
    if verbose: print('UI :',UL)
    if verbose: print('Hn :',Hn)
    if verbose: print('Ht :',Ht)
    if verbose: print('b  :',Mat)
    if verbose: print('dt :',dt)

    # element parameters
    NoElmtNodes = 4
    NoNodalDOF  = 1
    NoDim       = 2
    
    # initialize element vector /matrix
    r_e = np.zeros(NoElmtNodes*NoNodalDOF)
    k_e = np.zeros((NoElmtNodes*NoNodalDOF, NoElmtNodes*NoNodalDOF))

    # geometry and dofs
    XI = XL.reshape(-1,NoDim) # reformat to a matrix field: XI = [[X1x, X1y] ,[X2x, X2y] ,[X3x, X3y] ,[X4x, X4y]]
    TI = UL
    

    # Material Parameters
    alpha_q = Mat[0]


    # Provide integration points
    a = 1/np.sqrt(3)
    EGP = np.array([[-a,-a,1],[a,-a,1],[a,a,1],[-a,a,1]])
    NoInt = len(EGP)

    # Start integration Loop
    for GP in range(NoInt):
        if verbose: print('GP: ',GP)
        xi, eta, wgp  = EGP[GP]

        # compute current shape functions
        SH0 = SH0_Q1(xi, eta)

        # compute mapping
        Jed = np.einsum('Ii,Ij->ij', XI ,SH0[:,1:3])
        detJ = np.linalg.det(Jed)
        if (detJ <= 0): raise NameError("Error unphysical mapping detected.")
        if verbose: print('detJ: ',detJ)
        Jed_inv = np.linalg.inv(Jed)
        
        # map shape function derivative
        SHP = np.copy(SH0)
        SHP[:,1:3] = np.einsum('Ij,ji->Ii', SH0[:,1:3], Jed_inv)

        # compute gradient of temperature / heat flux
        gradT = np.einsum('Ii,I->i'  , SHP[:,1:3], TI)
        if verbose: print('gradT: ')
        if verbose: print(gradT)
        q = - alpha_q * gradT

        # export element vector
        for I in range(NoElmtNodes):

            # integrate nodal right hand side and export
            r_e[I*NoNodalDOF+0] += q.dot(SHP[I,1:3]) * wgp * detJ

            for J in range(NoElmtNodes):

                # compute nodal stiffness matrix
                k_e[I*NoNodalDOF+0, J*NoNodalDOF+0] += - alpha_q * SHP[J,1:3].dot( SHP[I,1:3])* wgp * detJ

    return r_e, k_e
```

Appendix: Shape functions for bi-linear quadirlateral elements:
```python
def SH0_Q1(xi, eta):
    '''
    SH0_Q1(xi, eta) -> SH0
    Return a two dimensional array containing shape functions and derivatives for Q1 element.
    Usage: SH0( NodeNumber, SHPIndex)  
    with SHPIndex = {
        0 -> shape function, 
        1 -> derived shape function w.r.t. xi, 
        2 -> derived shape function w.r.t. eta
        } 
    '''
    return 1/4 * np.array([
            [(1.0-xi)*(1.0-eta), -(1.0-eta),  -(1.0-xi)],
            [(1.0+xi)*(1.0-eta),  (1.0-eta),  -(1.0+xi)],
            [(1.0+xi)*(1.0+eta),  (1.0+eta),   (1.0+xi)],
            [(1.0-xi)*(1.0+eta), -(1.0+eta),   (1.0-xi)]
        ], dtype=np.float64)
```