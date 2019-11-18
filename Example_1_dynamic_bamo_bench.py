# TEST for FiniteElement in coupled problems
# for the dynamic terms including inertia and damping

import numpy as np
import matplotlib.pyplot as plt

from FEM_1D_SOLVER import FEM1D
import elements.Elmt_BaMo_BaEn_Coupled_1D as ELEMENT


# Create FEM Instance
FEM = FEM1D(ELEMENT)
FEM.Add_Mesh(9.0,1)
FEM.Add_Material([5,1.2,10,1,0,0],"All")
FEM.Add_EBC("x==0","U",0)
FEM.Add_EBC("x>-1","T",0)
FEM.Add_NBC("x==9",0,1)
FEM.Analysis()

# define a loading function
def load(time):
  lam = 0.0;
  if time <= 10:
    lam = (time/10)
  if time > 10:
    lam = 1.0
  if time > 60:
    lam = 0.0
  return lam

# Lets prepare a time loop, with recoding the time and displacement
rec_t = []
rec_u = []
rec_tu = []
nStep, time, dt = 100 ,0.0, 1.0
for step in range(nStep):
  time += dt
  FEM.NextStep(time,load(time))
  print( FEM.NewtonIteration() )
  print( FEM.NewtonIteration() )
  u = FEM.NodalDof("x==9","U")
  rec_t.append(time)
  rec_u.append(u)


plt.plot(rec_t,rec_u)
plt.xlabel('t')
plt.ylabel('u')
plt.show()
