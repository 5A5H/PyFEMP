# 2D example linear elasticity cooks membrane

import numpy as np
import matplotlib.pyplot as plt

import PyFEM
import PyFEM.elements.Elmt_BaMo_2D as ELEMENT

FEM = PyFEM.FEM_Simulation(ELEMENT)
XI, Elem = PyFEM.msh_rec([0.0, 0.0], [10.0, 10.0], [1, 1], type='Q1')
FEM.Add_Mesh(XI, Elem)
FEM.Add_Material([2100, 0.0], "All")
FEM.Add_EBC("x==0",  "UX", 0)
FEM.Add_EBC("y==0",  "UY", 0)
FEM.Add_EBC("x==10",  "UX", 1)
#FEM.Add_NBC("x==48", "UY", 10)
#FEM.Add_NBC("x==48 and (y==60 or y==44)" , "UY", 5)

FEM.Analysis()

FEM.NextStep(1.0, 1.0)
FEM.verbose = True

r,k,i = FEM.FormLinearSystem()
#print(i)
#FEM.CallElement(0)



#print( FEM.NewtonIteration() )
#print( FEM.NewtonIteration() )

