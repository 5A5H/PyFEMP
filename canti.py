# 2D example linear elasticity dynamic cantilever

import numpy as np
import matplotlib.pyplot as plt

import PyFEM
import PyFEM.elements.Elmt_BaMo_2D_Dy as ELEMENT

FEM = PyFEM.FEM_Simulation(ELEMENT)
n, sig = 3, 2
XI, Elem = PyFEM.msh_conv_quad([0.0, 0.0], [10.0, 0.0], [10.0, 1.0], [0.0, 1.0], [10*n, n], type='Q1')
FEM.Add_Mesh(XI, Elem)
FEM.Add_Material([2100, 0.3, 0.1, 0.0], "All")
FEM.Add_EBC("x==0",  "UX", 0)
FEM.Add_EBC("x==0",  "UY", 0)
FEM.Add_NBC("x==10", "UY", sig/n)
FEM.Add_NBC("x==10 and (y==0 or y==1)" , "UY", 1/2 * sig/n)

FEM.Analysis()

# define a loading function
def load(time):
  lam = 0.0
  if time <= 5:
    lam = (time/5)
  if time > 5:
    lam = 0.0
  return lam

# tempral discretization
nStep, time, dt = 100, 0.0, 0.1
animation = True

# record array
uy_vs_t = np.zeros((nStep+1,2))

for step in range(nStep):
  time += dt
  FEM.NextStep(0.01, 1.0)
  FEM.NextStep(time,load(time))
  print( FEM.NewtonIteration() )
  print( FEM.NewtonIteration() )
  uy = FEM.NodalDof("x==10 and y==1", "UY")
  uy_vs_t[step+1] = np.array([time, uy])
  
  if animation:
    if (step==0): fig, ax = plt.subplots(1,2, figsize=(10.0, 5.0))
    ax[0].cla()
    postplot = FEM.ShowMesh(ax[0], deformedmesh=True, PostName="SigMises")
    ax[0].set_xlim(0, 10.5)
    ax[0].set_ylim(-4.5, 4.5)
    
    ax[1].cla()
    ax[1].plot(uy_vs_t[:step+1, 0], uy_vs_t[:step+1, 1])
    ax[1].set_xlabel('t- time in s')
    ax[1].set_ylabel('uy- in m')
    ax[1].set_xlim(0,nStep*dt)
    ax[1].set_ylim(-4.5, 4.5)
    ax[1].set_aspect(9/9.5)
    plt.pause(0.001)
  

plt.show()




