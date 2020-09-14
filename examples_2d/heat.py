# 

import numpy as np
import matplotlib.pyplot as plt

import PyFEMP
import PyFEMP.elements.LE_T_Q1 as ELEMENT

FEM = PyFEMP.FEM_Simulation(ELEMENT)
XI, Elem = PyFEMP.msh_rec([0.0, 0.0], [10.0, 10.0], [40, 40])
FEM.Add_Mesh(XI, Elem)
FEM.Add_Material([2100, 0.3, 1.0], "All")
FEM.Add_EBC("x<1 and y<1",  "T" , 0)
FEM.Add_EBC("x>9 and y>9",  "T" , 0)
FEM.Add_EBC("x<1 and y>9",  "T" , 10)
FEM.Add_EBC("x>9 and y<1",  "T" , 10)
FEM.Add_EBC("x>4.5 and x<5.5 and y>4.5 and y<5.5",  "T" , 10)

FEM.Add_EBC("x==0",  "UX", 0)
FEM.Add_EBC("y==0",  "UY", 0)


FEM.Analysis()


FEM.NextStep(1.0, 1.0)

print( FEM.NewtonIteration() )
print( FEM.NewtonIteration() )

fig, ax = plt.subplots(1,1, figsize=(6.0, 8.0))
postplot = FEM.ShowMesh(ax, deformedmesh=True, PostName="T")
cbar = fig.colorbar(postplot)
cbar.ax.set_ylabel('Temperature $\theta$')
plt.show()


