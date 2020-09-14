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
FEM.Add_EBC("x==-2",  "UX", 0)
FEM.Add_EBC("y==-2",  "UY", 0)
FEM.Add_EBC("x== 2",  "UX", 0.5)


FEM.Analysis()

FEM.NextStep(1.0, 1.0)

print( FEM.NewtonIteration() )
print( FEM.NewtonIteration() )


fig, ax = plt.subplots(1,1, figsize=(8.0, 8.0))
postplot = FEM.ShowMesh(ax, deformedmesh=True, PostName="SigMises")
cbar = fig.colorbar(postplot)
cbar.ax.set_ylabel('von Mises stress $\sigma_{VM}$')
plt.show()


