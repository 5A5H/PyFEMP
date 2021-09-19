import numpy as np
import matplotlib.pyplot as plt

import PyFEMP
import T1_J2 as ELEMENT

FEM = PyFEMP.FEM_Simulation(ELEMENT)
n, sig = 10, 100
XI, Elem = PyFEMP.msh_conv_quad([0.0, 0.0], [48.0, 44.0], [48.0, 60.0], [0.0, 44.0], [2*n, n], type='T1')
FEM.Add_Mesh(XI, Elem)
FEM.Add_Material([2100, 0.3, 220.0, 25.0], "All")
FEM.Add_EBC("x==0",  "UX", 0)
FEM.Add_EBC("x==0",  "UY", 0)
FEM.Add_NBC("x==48", "UY", (sig*16)/n)
FEM.Add_NBC("x==48 and (y==60 or y==44)" , "UY", 1/2 * (sig*16)/n)

FEM.Analysis()

no_steps = 5
for step in range(no_steps):        
    FEM.NextStep((step+1), (step+1)*(1.0/no_steps))
    for i in range(6): residual = FEM.NewtonIteration()
    if (residual>1e-6):
        print("divergence in step: ", step, " with residual: ", residual)
        break

fig, axs = plt.subplots(1, 2, figsize=(16.0, 8.0))
postplot = FEM.ShowMesh(axs[0], deformedmesh=True, PostName="SigMises")
cbar = fig.colorbar(postplot, ax=axs[0])
cbar.ax.set_ylabel('von Mises stress $\sigma_{VM}$')
postplot2 = FEM.ShowMesh(axs[1], deformedmesh=True, PostName="a")
cbar2 = fig.colorbar(postplot2, ax=axs[1])
cbar2.ax.set_ylabel('eq. plastic arc length')
plt.show()