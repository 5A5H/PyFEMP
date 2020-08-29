# 2D example linear elasticity cooks membrane

import numpy as np
import matplotlib.pyplot as plt

import PyFEM
import PyFEM.elements.Elmt_BaMo_2D as ELEMENT

uu = []
nn = np.array([1, 2, 4, 6, 8, 10, 15, 20])
for n in nn:
    FEM = PyFEM.FEM_Simulation(ELEMENT)
    sig = 100
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

    uu.append([FEM.NoEquations, ux, uy])

uu = np.array(uu)
plt.plot(uu[:,0], -uu[:,1], label='ux')
plt.plot(uu[:,0],  uu[:,2], label='uy')
plt.title('Cook-Membrane Convergency study')
plt.xlabel('NoEquations')
plt.ylabel('displacement (upper right corner)')
plt.legend()
plt.show()