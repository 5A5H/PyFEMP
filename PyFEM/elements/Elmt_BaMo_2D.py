#from numpy import *
import numpy as np


def Elmt_Init():
    NoElementDim = 2
    NoElementNodes = 4
    NoElementHistory = 0
    ElementDofNames = ["UX", "UY"]
    ElementMaterialNames = ["E", "nu"]
    ElementPostNames = ["Sig"]
    return NoElementDim, NoElementNodes, ElementDofNames, NoElementHistory, ElementMaterialNames, ElementPostNames


def Elmt_KS(XI, UI, Hn, Ht, Mat, dt):
    '''
    '''
    verbose = False  # True;
    # element vector /matrix
    r_e = np.zeros(4*2)
    k_e = np.zeros((4*2, 4*2))
    # nodal coordinates
   
    return r_e, k_e


def Elmt_Post(XI, UI, Hn, Ht, Mat, dt, PostName):
    '''
    '''
    if PostName == "Sig":
        return 0, 0, 0, 0

    return X1, X2, 0, 0
