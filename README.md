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

### FEM_Simulation properties
Persistant Data provided by the FEM_Simulation object.
```
FEM_Simulation.NoElementDim                  # dimensions of this simulation
FEM_Simulation.NoElementNodes                # number of nodes for each element in the current simulation
FEM_Simulation.ElementDofNames               # vector with strings of names for nodal degrees of freedom
FEM_Simulation.NoElementHistory              # length of element history fields in the current simulation
FEM_Simulation.ElementMaterialNames          # vector with strings of material parameter names           
FEM_Simulation.ElementPostNames              # vector with strings of postprocessing names
FEM_Simulation.NoElementMaterial             # number of material parameters for each element
FEM_Simulation.NoNodeDofs                    # number of degrees of freedom per node in the current simulation

# general program variables
FEM_Simulation.verbose                        # verbose flag
FEM_Simulation.verbose_system                 # verbose flag
FEM_Simulation.state                          # current simulation state identifier

# general discretization variables
FEM_Simulation.time                          # current time
FEM_Simulation.dt                            # time increment gone frome last time
FEM_Simulation.step                          # current step
FEM_Simulation.lambda_load                   # global load multiplier
FEM_Simulation.NoElements                    # number of elements
FEM_Simulation.NoNodes                       # number of nodes
FEM_Simulation.NoDofs                        # number of degrees of freedom
FEM_Simulation.XI                            # nodal coordinates
FEM_Simulation.ELEM                          # element connectivity
FEM_Simulation.h_n                           # previous history field
FEM_Simulation.h_t                           # current history field

# initialize fields for boundary conditions
FEM_Simulation.NBC = []                       # python list to collect natural boundary conditions before analysis
FEM_Simulation.NBC_Indexes = 0                # vector of indexes to the external load vector where a nbc is present
FEM_Simulation.NBC_Values = 0                 # vector of values to be placed in the external load vector for each nbc index
FEM_Simulation.EBC = []                       # python list to collect essential boundary conditions before analysis
FEM_Simulation.EBC_Indexes = 0                # vector of indexes of constrained degrees of freedom
FEM_Simulation.EBC_Values = 0                 # vector of values for each constrained degree of freedom
FEM_Simulation.NoEquations = 0                # number of all unconstrained dofs

# element discretization parameter
FEM_Simulation.ElementMaterial = []           # list of material parameter
FEM_Simulation.h_n = 0                        # vector of element history field of t=t   (previous)
FEM_Simulation.h_t = 0                        # vector of element history field of t=t+1 (current)
FEM_Simulation.DI = 0                         # vector of degrees of freedom
FEM_Simulation.R_ext = 0                      # vector of external forces
```

### FEM_Simulation functions

**CallElementPost**
Calls the postprocessing routine Elmt_Post for element i with the current Simulation fields.
Returns first a list of all elment node indexes and next the vector with the requested PostName
data, one scalar for each node.
```
CallElementPost(self, i, PostName) -> elmt_nodes, r_post_e
```

**PostProcessing**
Returns a list of nodes and a list of 
one requested scalar for each of these nodes.
The requested scalar is computed from the element 
subroutine "Elmt_Post" by PostName.
By default, the PostName field is returned for all mesh nodes.
Optionally, specify arbitrary positions for which to evaluate the PostNames.
x -> matrix that contains the nodal positions
p -> vector of values for each node
```
PostProcessing("UX") -> x, p
PostProcessing("UX", [0.0, 0.0]) -> [0.0, 0.0], p
PostProcessing("UX", [[0.0, 0.0],...,[1.0, 0.0]]) -> [[0.0, 0.0],...,[1.0, 0.0]], p
```


## Using external meshes
External meshes can be used. Required is the input to the Add_Mesh command, as schown in the example ```plate.py```.
For PyFEMP it does not matter where the data come from,
once they are available in python they can be used.
**Attention**: The element connectivity requires indexing starting with **0**.
An example on exporting a T1 mesh from AceFEM/Mathematica is:
```
<< AceFEM`;
SMTInputData[];
SMTAddDomain[{"\[CapitalOmega]", {"ML:", "SE", "PE", "T1", "DF", "LE",
     "T1", "D", "Hooke"}, {"E *" -> 1000, "\[Nu] *" -> 0.3}}];
mesh = ToElementMesh[
   ImplicitRegion[x^2 + y^2 > 0.5, {x, y}], {{-1, 1}, {-1, 1}}, 
   "MeshOrder" -> 1, MaxCellMeasure -> 3, AccuracyGoal -> 1];
SMTAddMesh[mesh, "\[CapitalOmega]"];
SMTAnalysis[];
SMTShowMesh[]
XI = SMTNodeData[All, "X"];
Elmt = SMTElements[[;; , 3]];
Export["path/nodes.csv", XI];
Export["path/elmt.csv", Elmt - 1];
```
Notice the ```-1``` in the export for the element file, which translates to **0** indexing.

## Developer Notes
How to build the PyPi package from source, and upload/update it on PyPi to make it available via pip:

A version difference need to be specified in the ```setup.py``` file.
Currently, the files included in the package are detected automatically during build, as part of the standard python module structure.
The commands to re-build the package are:
```
python -m build
```
Now the packages are created.
To upload the package do:
```
python -m twine upload dist/*
```
and use the PyPi login as requested.

The whole process sometimes require deleting previous build folders or de-installation of previous PyFEMP version.

It is also possible to install the local repository after it has been build via
```
python -m pip install .
```
