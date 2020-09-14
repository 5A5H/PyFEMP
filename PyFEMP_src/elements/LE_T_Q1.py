import numpy as np


def Elmt_Init():
    NoElementDim = 2
    NoElementNodes = 4
    NoElementHistory = 0
    ElementDofNames = ["UX", "UY", "T"]
    ElementMaterialNames = ["E", "nu", "a"]
    ElementPostNames = ["UX", "UY", "T", "SigMises"]
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

def BmatVoigt2D(SHP):
    '''
    BmatVoigt(SHP) -> Bmat
    Returns a B-Matrix (as dim:3) for computing the strain vector in voigt notation.
    This B-Matrix assumes a 2D plane strain approximation.
    The is a shape function matrix with derived functions w.r.t. physical space.
    Input:
        SHP( NodeNumber, SHPIndex)  
        with SHPIndex = {
            0 -> shape function, 
            1 -> derived shape function w.r.t. x, 
            2 -> derived shape function w.r.t. y
            } 
    Output:
        Bmat(NodeNumber, i, j)
        for eps_i = B_Iij * u_Ij
        with 
        eps_i = [eps_11, eps_22, eps_33, 2*eps_12, 2*eps_23, 2*eps_13]
        u_Ij  = [[u_11, u_12] ... [u_n1, u_n2]]
    '''
    return np.array([ 
        [
            [N[1], 0   ],
            [0   , N[2]],
            [0   , 0   ],
            [N[2], N[1]],
            [0   , 0   ],
            [0   , 0   ]
        ]
        for N in SHP 
     ], dtype=np.float64)

def HookeMatVoigt(lam, mue):
    '''
    HookeMatVoigt(lam, mue) -> Cmat
    Returns the constitutive Voigt MATRIX(6,6) for a Hooke material law.
    The input are the elastic Lame constants.
    sig_i = Cmat_ij * eps_j
    with
    sig_i = [sig_11, sig_22, sig_33, sig_12, sig_23, sig_13]
    eps_i = [eps_11, eps_22, eps_33, 2*eps_12, 2*eps_23, 2*eps_13]
    '''
    return np.array([
        [lam + 2* mue, lam         , lam         , 0  , 0  , 0  ],
        [lam         , lam + 2* mue, lam         ,0   , 0  , 0  ],
        [lam         , lam         , lam + 2* mue, 0  , 0  , 0  ],
        [0           , 0           , 0           , mue, 0  , 0  ],
        [0           , 0           , 0           , 0  , mue, 0  ],
        [0           , 0           , 0           , 0  , 0  , mue]
    ], dtype=np.float64)        
    

def Elmt_KS(XL, UL, Hn, Ht, Mat, dt):
    '''
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
    NoNodalDOF  = 3
    NoDim       = 2
    
    # initialize element vector /matrix
    r_e = np.zeros(NoElmtNodes*NoNodalDOF)
    k_e = np.zeros((NoElmtNodes*NoNodalDOF, NoElmtNodes*NoNodalDOF))

    # geometry and dofs
    XI = XL.reshape(-1,NoDim)
    uxI = UL[0::NoNodalDOF]
    uyI = UL[1::NoNodalDOF]
    TI  = UL[2::NoNodalDOF]
    uI = np.array([uxI, uyI]).T
    

    # Material Parameters
    Emod, nu = Mat[0], Mat[1]
    lam , mue = (Emod*nu)/((1.0+nu)*(1.0-2.0*nu)), Emod/(2.0*(1.0+nu))

    # constitutive matrix (hooke)
    Cmat = HookeMatVoigt(lam, mue)

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
        Bmat = BmatVoigt2D(SHP)
        # compute strains / stresses
        eps = np.einsum('Iij,Ij->i', Bmat, uI)
        if verbose: print('eps: ')
        if verbose: print(np.array([eps[0], eps[3]/2, eps[4]/2, eps[3]/2, eps[1], eps[5]/2, eps[4]/2, eps[5]/2, eps[2]]))
        sig = np.einsum('ij,j->i'  , Cmat, eps)
        if verbose: print('sig: ')
        if verbose: print(np.array([[sig[0], sig[3], sig[4]],[sig[3], sig[1], sig[5]],[sig[4], sig[5], sig[2]]]))
        # compute gradient of temperature / heat flux
        gradT = np.einsum('Ii,I->i'  , SHP[:,1:3], TI)
        if verbose: print('gradT: ')
        if verbose: print(gradT)
        q = - a * gradT

        # export right hand side | this element has 4 nodes with 2 dofs each
        for I in range(NoElmtNodes):
            
            # compute nodal right hand side
            nodal_rhs_vec    = np.einsum('i,ij->j',sig, Bmat[I])

            # integrate nodal right hand side and export
            r_e[I*NoNodalDOF+0] += nodal_rhs_vec[0] * wgp * detJ
            r_e[I*NoNodalDOF+1] += nodal_rhs_vec[1] * wgp * detJ
            r_e[I*NoNodalDOF+2] += q.dot(SHP[I,1:3]) * wgp * detJ

            for J in range(NoElmtNodes):

                # compute nodal stiffness matrix
                nodal_stiffness = np.einsum('ki,ko,oj->ij', Bmat[I], Cmat, Bmat[J])
                k_e[I*NoNodalDOF+0, J*NoNodalDOF+0] += nodal_stiffness[0,0] * wgp * detJ
                k_e[I*NoNodalDOF+0, J*NoNodalDOF+1] += nodal_stiffness[0,1] * wgp * detJ
                k_e[I*NoNodalDOF+1, J*NoNodalDOF+0] += nodal_stiffness[1,0] * wgp * detJ
                k_e[I*NoNodalDOF+1, J*NoNodalDOF+1] += nodal_stiffness[1,1] * wgp * detJ

                k_e[I*NoNodalDOF+2, J*NoNodalDOF+2] += - a * SHP[J,1:3].dot( SHP[I,1:3])* wgp * detJ


    return r_e, k_e


def Elmt_Post(XL, UL, Hn, Ht, Mat, dt, PostName):
    '''
    '''
    ## NEW post strategy is to return indexes
    # and contribuzions
    
    # element parameters
    NoElmtNodes = 4
    NoNodalDOF  = 3
    NoDim       = 2

    # 4 nodes
    r_post = np.zeros(NoElmtNodes)

    # geometry and dofs
    XI = XL.reshape(-1,NoDim)
    uxI = UL[0::NoNodalDOF]
    uyI = UL[1::NoNodalDOF]
    TI  = UL[2::NoNodalDOF]
    uI = np.array([uxI, uyI]).T
    

    # Material Parameters
    Emod, nu = Mat[0], Mat[1]
    lam , mue = (Emod*nu)/((1.0+nu)*(1.0-2.0*nu)), Emod/(2.0*(1.0+nu))

    # constitutive matrix (hooke)
    Cmat = HookeMatVoigt(lam, mue)

    # Provide integration points
    a = 1/np.sqrt(3)
    EGP = np.array([[-a,-a,1],[a,-a,1],[a,a,1],[-a,a,1]])
    NoInt = len(EGP)

    # Start integration Loop
    for GP in range(NoInt):
        xi, eta, wgp  = EGP[GP]

        # compute current shape functions
        SH0 = SH0_Q1(xi, eta)

        # compute mapping
        Jed = np.einsum('Ii,Ij->ij', XI ,SH0[:,1:3])
        detJ = np.linalg.det(Jed)
        if (detJ <= 0): raise NameError("Error unphysical mapping detected.")
        Jed_inv = np.linalg.inv(Jed)
        
        # map shape function derivative
        SHP = np.copy(SH0)
        SHP[:,1:3] = np.einsum('Ij,ji->Ii', SH0[:,1:3], Jed_inv)
        Bmat = BmatVoigt2D(SHP)

        # compute strains / stresses
        eps = np.einsum('Iij,Ij->i', Bmat, uI)
        sig = np.einsum('ij,j->i'  , Cmat, eps)

        sig_vm = np.sqrt( \
            sig[0]**2 + sig[1]**2 + sig[2]**2 \
            -sig[0]*sig[1] -sig[0]*sig[2] -sig[1]*sig[2] \
            + 3* (sig[3]**2 + sig[4]**2 + sig[5]**2) \
            )
        
        if (PostName=="SigMises"):
            r_post += sig_vm * SHP[:,0]

    if (PostName=="UX"):
        r_post = np.array([UL[0], UL[3], UL[6], UL[9]])
        return r_post
    elif (PostName=="UY"):
        r_post = np.array([UL[1], UL[4], UL[7], UL[10]])
        return r_post
    elif (PostName=="T"):
        r_post = np.array([UL[2], UL[5], UL[8], UL[11]])
        return r_post
    elif (PostName=="SigMises"):
        return r_post
    else:
        print("Waring: PostName "+PostName+" not defined!")
        return np.array([0.0, 0.0, 0.0, 0.0])
