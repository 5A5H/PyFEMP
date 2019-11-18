# TEST for FiniteElement in coupled problems 
# for the dynamic terms including inertia and damping 

import numpy as np
import matplotlib.pyplot as plt
from FEM_1D_SOLVER import FEM1D
import Elmt_BaMo_BaEn_Coupled_1D as ELEMENT

# Create FEM Instance
FEM = FEM1D(ELEMENT)
FEM.verbose_system = False
FEM.Add_Mesh(100.0,100)
FEM.Add_Material([100,1,1,100,0,0],"All")
FEM.Mod_Material([0,1,1,100,0,0],60)
FEM.Add_EBC("x==0","U",0)
FEM.Add_EBC("x>-1","T",0)
FEM.Add_NBC("x==100",0,1)
FEM.Analysis()

def PlotDomain():
  x = FEM.XI
  u = FEM.DI[0::2]
  plt.clf()
  plt.plot(x,u)
  plt.xlabel('x')
  plt.ylabel('u')
  plt.xlim(0,100)
  plt.ylim(-0.2,0.2)
  plt.title('t=%f4'%FEM.time)
  plt.show()  

# define a loading function
def load(time):
  lam = 0.0;
  if time <= 5:
    lam = (time/5)
  if time > 5:
    lam = 1-((time-5)/5)
  if time > 10:
    lam = 0.0
  return lam

# Lets prepare a time loop, with recoding the time and displacement
rec_t = []
rec_u = []
#rec_tu = []
nStep, time, dt = 200 ,0.0, 1.0
for step in range(nStep):
  time += dt
  FEM.NextStep(time,load(time))
  FEM.NewtonIteration()
  FEM.NewtonIteration()
  PlotDomain()
#  u = FEM.NodalDof("x==9","U")
  rec_t.append(time)
  rec_u.append(load(time))


plt.plot(rec_t,rec_u)
#plt.xlabel('t')
#plt.ylabel('u')
plt.show()
