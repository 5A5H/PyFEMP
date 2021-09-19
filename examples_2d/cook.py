# 2D example linear elasticity cooks membrane

import numpy as np
import matplotlib.pyplot as plt

import PyFEMP
import PyFEMP.elements.Elmt_BaMo_2D as ELEMENT

FEM = PyFEMP.FEM_Simulation(ELEMENT)
n, sig = 2, 100
XI, Elem = PyFEMP.msh_conv_quad([0.0, 0.0], [48.0, 44.0], [48.0, 60.0], [0.0, 44.0], [2*n, n], type='Q1')
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

fig, axs = plt.subplots(1, 2, figsize=(12.0, 8.0))
postplot = FEM.ShowMesh(axs[0], boundaryconditions=True)
axs[0].legend()
postplot = FEM.ShowMesh(axs[1], deformedmesh=True, PostName="SigMises")
cbar = fig.colorbar(postplot)
cbar.ax.set_ylabel('von Mises stress $\sigma_{VM}$')
plt.show()




