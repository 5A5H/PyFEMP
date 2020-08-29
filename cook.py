# 2D example linear elasticity cooks membrane

import numpy as np
import matplotlib.pyplot as plt

import PyFEM
import PyFEM.elements.Elmt_BaMo_2D as ELEMENT

FEM = PyFEM.FEM_Simulation(ELEMENT)
n, sig = 4, 100
XI, Elem = PyFEM.msh_conv_quad([0.0, 0.0], [48.0, 44.0], [48.0, 60.0], [0.0, 44.0], [2*n, n], type='Q1')
FEM.Add_Mesh(XI, Elem)
FEM.Add_Material([2100, 0.3], "All")
FEM.Add_EBC("x==0",  "UX", 0)
FEM.Add_EBC("x==0",  "UY", 0)
FEM.Add_NBC("x==48", "UY", (sig*16)/n)
FEM.Add_NBC("x==48 and (y==60 or y==44)" , "UY", 1/2 * (sig*16)/n)

FEM.Analysis()

FEM.NextStep(1.0, 1.0)
print( FEM.NewtonIteration() )
print( FEM.NewtonIteration() )

ux = FEM.NodalDof("x==48 and y==60","UX")
uy = FEM.NodalDof("x==48 and y==60","UY")
print('ux :',ux, 'uy :',uy)

fig, ax = plt.subplots(1,1, figsize=(6.0, 6.0))
#FEM.ShowMesh(ax)
postplot = FEM.ShowMesh(ax, deformedmesh=True, PostName="SigMises")
fig.colorbar(postplot)

plt.show()




