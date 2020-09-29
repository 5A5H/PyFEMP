import numpy as np


def Elmt_Init():
    '''
    The initialization function of a PyFEMP element.
    '''
    NoElementDim = 2                    # Number of dimensions
    NoElementNodes = 4                  # Number of nodes on the element (here: 4-Node Quadirlateral)
    NoElementHistory = 0                # Number of history variables
    ElementDofNames = ["T"]             # Element degrees of freedom specified by a string for each scalar
    ElementMaterialNames = ["alpha_q"]  # Element material parameters specified by a string for each scalar
    ElementPostNames = ["T"]            # Element postprocessing parameters
    return NoElementDim, NoElementNodes, ElementDofNames, NoElementHistory, ElementMaterialNames, ElementPostNames      
    
def SH0_Q1(xi, eta):
    '''
    SH0_Q1(xi, eta) -> SH0
    Return a two dimensional array containing shape functions and derivatives for Q1 element.
    Usage: SH0( NodeNumber, SHPIndex)  
    with SHPIndex = {
        0 -> shape function, 
        1 -> derived shape function w.r.t. xi, 
        2 -> derived shape function w.r.t. eta
        } 
    '''
    return 1/4 * np.array([
            [(1.0-xi)*(1.0-eta), -(1.0-eta),  -(1.0-xi)],
            [(1.0+xi)*(1.0-eta),  (1.0-eta),  -(1.0+xi)],
            [(1.0+xi)*(1.0+eta),  (1.0+eta),   (1.0+xi)],
            [(1.0-xi)*(1.0+eta), -(1.0+eta),   (1.0-xi)]
        ], dtype=np.float64)

def Elmt_KS(XL, UL, Hn, Ht, Mat, dt):
    '''
    This function returns element vector and matrix.
    Generally the element matrix is the straight derivative
    of the element vector with respect to the nodal degrees of freedoms,
    in the order of the UL field.
    Both need to be returned as numpy arrays, with dimensions: 
    element vector r_e : number of nodal degrees of freedom times number of nodes
    element matrix k_e : NoNodes * NoNodalDOF ,NoNodes * NoNodalDOF
    '''
    verbose = False  # True;
    if verbose: print('XI :',XL)
    if verbose: print('UI :',UL)
    if verbose: print('Hn :',Hn)
    if verbose: print('Ht :',Ht)
    if verbose: print('b  :',Mat)
    if verbose: print('dt :',dt)

    # element parameters
    NoElmtNodes = 4
    NoNodalDOF  = 1
    NoDim       = 2
    
    # initialize element vector /matrix
    r_e = np.zeros(NoElmtNodes*NoNodalDOF)
    k_e = np.zeros((NoElmtNodes*NoNodalDOF, NoElmtNodes*NoNodalDOF))

    # geometry and dofs
    XI = XL.reshape(-1,NoDim) # reformat to a matrix field: XI = [[X1x, X1y] ,[X2x, X2y] ,[X3x, X3y] ,[X4x, X4y]]
    TI = UL
    

    # Material Parameters
    alpha_q = Mat[0]


    # Provide integration points
    a = 1/np.sqrt(3)
    EGP = np.array([[-a,-a,1],[a,-a,1],[a,a,1],[-a,a,1]])
    NoInt = len(EGP)

    # Start integration Loop
    for GP in range(NoInt):
        if verbose: print('GP: ',GP)
        xi, eta, wgp  = EGP[GP]

        # compute current shape functions
        SH0 = SH0_Q1(xi, eta)

        # compute mapping
        Jed = np.einsum('Ii,Ij->ij', XI ,SH0[:,1:3])
        detJ = np.linalg.det(Jed)
        if (detJ <= 0): raise NameError("Error unphysical mapping detected.")
        if verbose: print('detJ: ',detJ)
        Jed_inv = np.linalg.inv(Jed)
        
        # map shape function derivative
        SHP = np.copy(SH0)
        SHP[:,1:3] = np.einsum('Ij,ji->Ii', SH0[:,1:3], Jed_inv)

        # compute gradient of temperature / heat flux
        gradT = np.einsum('Ii,I->i'  , SHP[:,1:3], TI)
        if verbose: print('gradT: ')
        if verbose: print(gradT)
        q = - alpha_q * gradT

        # export element vector
        for I in range(NoElmtNodes):

            # integrate nodal right hand side and export
            r_e[I*NoNodalDOF+0] += q.dot(SHP[I,1:3]) * wgp * detJ

            for J in range(NoElmtNodes):

                # compute nodal stiffness matrix
                k_e[I*NoNodalDOF+0, J*NoNodalDOF+0] += - alpha_q * SHP[J,1:3].dot( SHP[I,1:3])* wgp * detJ

    return r_e, k_e


def Elmt_Post(XL, UL, Hn, Ht, Mat, dt, PostName):
    '''
    '''
    if (PostName=="T"):
        r_post = np.array([UL[0], UL[1], UL[2], UL[3]])
        return r_post
    else:
        print("Waring: PostName "+PostName+" not defined!")
        return np.array([0.0, 0.0, 0.0, 0.0])
