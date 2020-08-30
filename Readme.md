# **PyFEMP** **Py**thon **F**inite **E**lement **P**rogram
 PyFEMP (Python Finite Element Program) is a simple Finite Element program written in python. Its focus is on simplicity **not** on performance.

![Canti](PyFEMP_src/assets/canti.png?raw=true "Dynamic analysis of a cantilever")

 It should be easy to use, to understand and as portable as possible to be used in teaching. We aim to void overhead w.r.t. environmental setup (compiler, libraries, e.t.c. ...) and dealing with complex structures, to focus on the essense of the FEM.

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
Once downloaded this project you just need to execute the setup script and should be ready to run:
```
python setup.py
```
>> Remark: For **updating** PyFEMP just download the new project and run the setup script again. It always overwrites!

More information will follow in the documentation.

# Usage
Once succesfully installed you can start by surveying and running the examples, e.g.
```
python cook.py
```

 ![Canti](PyFEMP_src/assets/cook.png?raw=true "Dynamic analysis of a cantilever")