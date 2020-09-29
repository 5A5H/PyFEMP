# **PyFEMP** **Py**thon **F**inite **E**lement **P**rogram
 PyFEMP (Python Finite Element Program) is a simple Finite Element program written in python. Its focus is on simplicity **not** on performance.

![Canti](PyFEMP/assets/canti.png?raw=true "Dynamic analysis of a cantilever")

 It should be easy to use, to understand and as portable as possible to be used in teaching. We aim to void overhead w.r.t. environmental setup (compiler, libraries, e.t.c. ...) or dealing with complex structures, to focus on the essense of the FEM.

 Therefore PyFEMP is written completely in python with the only required modules are *numpy* and *matplotlib*. Furthermore the program provides lsess than 30 commands, including postprocessing and meshing. Python code is really easy to read and write. This holds especially when performance is not critical. In PyFEMP we emprace python loops (easy to read/write but weak in performance), leading to acceptable runtime (<20s) up to ~2000 Elements.

 # What it can do:
 * static/dynamic FE analysis
 * 1D, 2D analysis (3D is techniqually possible but visualisation is not supported via matplotlib) 
 * implementing user specific elements
 * arbitrary number of nodal degrees of freedom e.g. for coupled FEM
 * nonlinear analysisi via implemented Newton-Raphson
 * visualsiation of element specific postprocessing fields
 * access to all variables e.g. element matrices during the simulation or the linear equation system
 * generation of simple regular meshes in 1D/2D using triangl- and quadrilateral elements (external meshes can be used)

 # Install

 ## Requirements
 Using PyFEMP requires a python (3.x) installation with *numpy* and *matplotlib*. 

 We recomend installing *Anaconda*. Its a free python distribution that comes with poth *numpy* and *matplotlib* installed by default and is available for all major computing platforms. 
 This should result in the **works out of the box** experience:
 >> https://www.anaconda.com/products/individual

## Installation
The package is available to be installed via pip:
```
pip install --upgrade PyFEMP
```

# Usage
Once succesfully installed you can start by surveying and running the examples, e.g.
```
python cook.py
```

 ![Canti](PyFEMP/assets/cook.png?raw=true "Dynamic analysis of a cantilever")

## Function Reference
* functions of `PyFEMP`:
   + `PyFEMP.FEM_Simulation(ELEMENT)`
      Used to start a FEM Simlation by returnig the simulation object.
      The input is an element, providing the nessercary functions.

   + `PyFEMP.msh_line(X0, X1, N, type='U1') -> XI, ELEM`
      Used to generate a regular mesh in 1D. Returns a list of nodal coordinatex and 
      a matrix of element connectivity.

   + `PyFEMP.msh_rec(X0, X1, N, type='Q1') -> XI, ELEM`
      Used to generate a regular mesh in 2D for a rectangle, specified
      by the lower left corner `X0` and upper right `X1`.
      Returns a list of nodal coordinatex and a matrix of element connectivity.

   + `PyFEMP.msh_conv_quad(X1, X2, X3, X4, N, type='Q1') -> XI, ELEM`
      Same as `msh_rec` but for an arbitrary convex quadrilateral, specified
      by vertices `X1, X2, X3, X4`.

## The `PyFEMP.FEM_Simulation` object
The `PyFEMP.FEM_Simulation` object represents your FEM Simulation. It provides the
methods to e.g. introduce the mesh and boundary conditions i.e the discretized boundary 
value problem, but also to perform solution procedures on it.

### properties

