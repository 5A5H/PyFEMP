# reading meshes from csv

import os
import numpy as np
import matplotlib.pyplot as plt

import PyFEMP
import PyFEMP.elements.LE_T1 as ELEMENT

FEM = PyFEMP.FEM_Simulation(ELEMENT)
# Here we chose one of three pais of node and element files
# shipped with PyFEMP
Mesh = [
    ["plate_with_hole_nodes.csv"           , "plate_with_hole_elmts.csv"],
    ["plate_with_hole_nodes_fine.csv"      , "plate_with_hole_elmts_fine.csv"],
    ["plate_with_hole_nodes_super_fine.csv", "plate_with_hole_elmts_super_fine.csv"]
][2]

CSV_NodeFile = os.path.join(PyFEMP.assets, "plate_with_hole", Mesh[0])
CSV_ElmtFile = os.path.join(PyFEMP.assets, "plate_with_hole", Mesh[1])
XI   = np.genfromtxt(CSV_NodeFile, delimiter=',')
Elem = np.genfromtxt(CSV_ElmtFile, delimiter=',')
FEM.Add_Mesh(XI, Elem)
FEM.Add_Material([2100, 0.3], "All")
# external meshes might require more loose conditionals 
# as rounding errors during io might occur
FEM.Add_EBC("x<=-1", "UX", 0)
FEM.Add_EBC("x>=1",  "UX", 0.5)
FEM.Add_EBC("y==-1 and x==-1",  "UY", 0)


FEM.Analysis()

FEM.NextStep(1.0, 1.0)

print( FEM.NewtonIteration() )
print( FEM.NewtonIteration() )

# examples for 2d postprocessing:
# requesting post-data for all mesh-nodes
allnodes, postforallnodes = FEM.PostProcessing("SigMises")

# requesting post-data for a specific point (reference configuration only) (PostName must be implemented in the element)
requestednodes, postforrequestednodes = FEM.PostProcessing("SigMises", [0.0, 0.0])
print('Stress at ', requestednodes, '  ->  ', postforrequestednodes)
requestednodes, postforrequestednodes = FEM.PostProcessing("SigMises", [-0.9, 0.12])
print('Stress at ', requestednodes, '  ->  ', postforrequestednodes)

# requesting post data for a set of points and plot them
x_diag = np.linspace(-1, 1, 30)
y_diag = np.linspace(-1, 1, 30).T
requestednodes = np.array([x_diag, y_diag]).T
requestednodes, postforrequestednodes = FEM.PostProcessing("SigMises", requestednodes)
requestednodes_projected = np.linalg.norm(requestednodes-[-1,-1],axis=1)

fig, axs = plt.subplots(1, 2, figsize=(16.0, 8.0))
postplot = FEM.ShowMesh(axs[0], deformedmesh=True, PostName="SigMises")
cbar = fig.colorbar(postplot)
cbar.ax.set_ylabel('von Mises stress $\sigma_{VM}$')
axs[1].plot(requestednodes_projected, postforrequestednodes)
axs[1].fill_between(requestednodes_projected, postforrequestednodes)
axs[1].set_xlabel('x- distance from [0.0, 0.0]')
axs[1].set_ylabel('von Mises stress $\sigma_{VM}$')
plt.show()