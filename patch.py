# 2D example tensile test

import numpy as np
import matplotlib.pyplot as plt

import PyFEM
import PyFEM.elements.Elmt_BaMo_2D as ELEMENT

FEM = PyFEM.FEM_Simulation(ELEMENT)
n = 4
XI, Elem = PyFEM.msh_rec([0.0, 0.0], [10.0, 10.0], [n, n], type='Q1')
FEM.Add_Mesh(XI, Elem)
FEM.Add_Material([2100, 0.3], "All")
FEM.Add_EBC("x==0",  "UX", 0)
FEM.Add_EBC("y==0",  "UY", 0)
FEM.Add_EBC("x==10",  "UX", 1)

FEM.Analysis()

FEM.NextStep(1.0, 1.0)

print( FEM.NewtonIteration() )
print( FEM.NewtonIteration() )

ux = FEM.NodalDof("x==10 and y==10", "UX")
uy = FEM.NodalDof("x==10 and y==10", "UY")
print('ux :',ux, 'uy :',uy)

fig, ax = plt.subplots(1,1, figsize=(8.0, 8.0))
postplot = FEM.ShowMesh(ax, ec='b', label='reference config.')
postplot = FEM.ShowMesh(ax, deformedmesh=True, ec='r', label='current config.')
ax.legend()
plt.show()