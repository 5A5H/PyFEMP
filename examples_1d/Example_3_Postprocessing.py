# TEST for FiniteElement in coupled problems
# postprocessing of a truss under gravity
# seperately solving laplace -> temperature in the middle and on one side set

import matplotlib.pyplot as plt

import PyFEM
import PyFEM.elements.Elmt_BaMo_BaEn_Coupled_1D as ELEMENT


# Create FEM Instance
FEM = PyFEM.FEM_Simulation(ELEMENT)
FEM.Add_1DMesh(10.0,20)
FEM.Add_Material([5,1.2,0,0,1,0],"All")
FEM.Add_EBC("x==0","U",0)
FEM.Add_EBC("x==10","T",0)
FEM.Add_EBC("x==5","T",3)
for node in range(FEM.NoNodes):
  FEM.Add_NBC(node, "U", 0.1)

FEM.Analysis()

# Analysis of a static problem
FEM.NextStep(1,1)
FEM.NewtonIteration()
FEM.NewtonIteration()

#Plot Accelerations,Stresses over 1D domain
plt.figure(1,figsize=[20,8])

XI = FEM.XI
UI = FEM.DI[0::2]
plt.subplot(221)
plt.plot(XI,UI)
plt.xlabel('x')
plt.ylabel('$u$')

XI = FEM.XI
TI = FEM.DI[1::2]
plt.subplot(222)
plt.plot(XI,TI)
plt.xlabel('x')
plt.ylabel('$T$')

XI, SigI = FEM.PostProcessing("Sig")
plt.subplot(223)
plt.plot(XI,SigI)
plt.xlabel('x')
plt.ylabel('$\sigma$')

XI, A = FEM.PostProcessing("q")
plt.subplot(224)
plt.plot(XI,A)
plt.xlabel('x')
plt.ylabel('q')

plt.show()