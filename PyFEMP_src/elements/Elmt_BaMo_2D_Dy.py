import numpy as np


def Elmt_Init():
    NoElementDim = 2
    NoElementNodes = 4
    NoElementHistory = 4*2*3
    ElementDofNames = ["UX", "UY"]
    ElementMaterialNames = ["E", "nu", "rho", "d"]
    ElementPostNames = ["UX", "UY", "SigMises"]
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
        [lam         , lam + 2* mue, lam         , 0  , 0  , 0  ],
        [lam         , lam         , lam + 2* mue, 0  , 0  , 0  ],
        [0           , 0           , 0           , mue, 0  , 0  ],
        [0           , 0           , 0           , 0  , mue, 0  ],
        [0           , 0           , 0           , 0  , 0  , mue]
    ], dtype=np.float64)

def DampMatVoigt(damp):
    '''
    DampMatVoigt(damp) -> Dmat
    Returns the constitutive Voigt MATRIX(6,6) for a isotropic damping law.
    The input is a single damping parameter.
    sig_i_damp = Dmat_ij * deps_j
    with
    sig_i  = [sig_11 , sig_22 , sig_33 , sig_12   , sig_23   , sig_13   ]
    deps_i = [deps_11, deps_22, deps_33, 2*deps_12, 2*deps_23, 2*deps_13]
    '''
    return np.array([
        [damp  , 0     , 0     , 0     , 0     , 0     ],
        [0     , damp  , 0     , 0     , 0     , 0     ],
        [0     , 0     , damp  , 0     , 0     , 0     ],
        [0     , 0     , 0     , damp/2, 0     , 0     ],
        [0     , 0     , 0     , 0     , damp/2, 0     ],
        [0     , 0     , 0     , 0     , 0     , damp/2]
    ], dtype=np.float64)          

def Newmark_V(U, U_n, V_n, A_n, gamma, beta, dt):
    return gamma/(beta*dt) * (U-U_n) + (1.0-gamma/beta) * V_n + dt * (1.0-gamma/(2.0*beta)) * A_n


def Newmark_A(U, U_n, V_n, A_n, gamma, beta, dt):
    return 1.0/(beta*dt*dt) * (U-U_n - dt*V_n - dt*dt*(1.0/2.0 - beta)*A_n)
    

def Elmt_KS(XL, UL, Hn, Ht, Mat, dt):
    '''
    '''
    verbose = False  # True;
    if verbose: print('XI :',XL)    # XL  = [x11, x12, x21, x22, x31, x32, x41, x42]
    if verbose: print('UI :',UL)    # UL  = [u11, u12, u21, u22, u31, u32, u41, u42]
    if verbose: print('Hn :',Hn)    # Hn  = [a11, a12, v11, v12, ..., a41, a42, v41, v42]
    if verbose: print('Ht :',Ht)    # Ht  = [a11, a12, v11, v12, ..., a41, a42, v41, v42]
    if verbose: print('b  :',Mat)   # Mat = [Emod, nue, rho, damp]
    if verbose: print('dt :',dt)    # dt  = dt
    
    # initialize element vector /matrix
    r_e = np.zeros(4*2)
    k_e = np.zeros((4*2, 4*2))

    # geometry and dofs
    XI = np.array([[XL[0], XL[1]], [XL[2], XL[3]], [XL[ 4], XL[ 5]], [XL[ 6], XL[ 7]]], dtype=np.float64)
    uI = np.array([[UL[0], UL[1]], [UL[2], UL[3]], [UL[ 4], UL[ 5]], [UL[ 6], UL[ 7]]], dtype=np.float64)
    
    # read in history
    aIn = np.array([[Hn[ 0], Hn[ 1]], [Hn[ 6], Hn[ 7]], [Hn[12], Hn[13]], [Hn[18], Hn[19]]], dtype=np.float64)
    vIn = np.array([[Hn[ 2], Hn[ 3]], [Hn[ 8], Hn[ 9]], [Hn[14], Hn[15]], [Hn[20], Hn[21]]], dtype=np.float64)
    uIn = np.array([[Hn[ 4], Hn[ 5]], [Hn[10], Hn[11]], [Hn[16], Hn[17]], [Hn[22], Hn[23]]], dtype=np.float64)

    # newmark
    gamma = 1.0/2.0
    beta  = 1.0/4.0
    aI    = Newmark_A(uI, uIn, vIn, aIn, gamma, beta, dt)
    vI    = Newmark_V(uI, uIn, vIn, aIn, gamma, beta, dt)

    # export history
    [[Ht[ 0], Ht[ 1]], [Ht[ 6], Ht[ 7]], [Ht[12], Ht[13]], [Ht[18], Ht[19]]] = aI
    [[Ht[ 2], Ht[ 3]], [Ht[ 8], Ht[ 9]], [Ht[14], Ht[15]], [Ht[20], Ht[21]]] = vI
    [[Ht[ 4], Ht[ 5]], [Ht[10], Ht[11]], [Ht[16], Ht[17]], [Ht[22], Ht[23]]] = uI
    

    # Material Parameters
    Emod, nu, rho , d = Mat[0], Mat[1], Mat[2], Mat[3]
    lam , mue = (Emod*nu)/((1.0+nu)*(1.0-2.0*nu)), Emod/(2.0*(1.0+nu))

    # constitutive matrix (hooke)
    Cmat = HookeMatVoigt(lam, mue)

    # damping matrix (isotropic)
    Dmat = DampMatVoigt(d)

    # provide integration points
    aa = 1/np.sqrt(3)
    EGP = np.array([[-aa, -aa, 1],[aa, -aa, 1],[aa, aa, 1],[-aa, aa, 1]])
    NoInt = len(EGP)

    # start integration Loop
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

        deps = np.einsum('Iij,Ij->i', Bmat, vI)
        if verbose: print('deps: ')
        if verbose: print(np.array([deps[0], deps[3]/2, deps[4]/2, deps[3]/2, deps[1], deps[5]/2, deps[4]/2, deps[5]/2, deps[2]]))
        depsdu = Bmat * (gamma/(beta*dt))

        sig = np.einsum('ij,j->i'  , Cmat, eps) \
            + np.einsum('ij,j->i'  , Dmat, deps)
        if verbose: print('sig: ')
        if verbose: print(np.array([[sig[0], sig[3], sig[4]],[sig[3], sig[1], sig[5]],[sig[4], sig[5], sig[2]]]))

        # compute acelleration
        a = np.einsum('I,Ii->i',SHP[:,0], aI)
        dadu = (1.0/(beta*dt**2)) * np.eye(2)

        # export right hand side | this element has 4 nodes with 2 dofs each
        for I in range(4):
            
            # compute nodal right hand side
            nodal_rhs_vec    = np.einsum('i,ij->j',sig, Bmat[I]) + a * SHP[I,0] * rho

            # integrate nodal right hand side and export
            r_e[I*2+0] += nodal_rhs_vec[0] * wgp * detJ
            r_e[I*2+1] += nodal_rhs_vec[1] * wgp * detJ

            for J in range(4):

                # compute nodal stiffness matrix
                nodal_stiffness = np.einsum('ki,ko,oj->ij', Bmat[I], Cmat, Bmat[J]) \
                                + SHP[I,0] * SHP[J,0] * rho * dadu \
                                + np.einsum('ki,ko,oj->ij', depsdu[I], Dmat, Bmat[J])

                # integrate nodal stiffness matrix and export
                k_e[I*2+0, J*2+0] += nodal_stiffness[0,0] * wgp * detJ
                k_e[I*2+0, J*2+1] += nodal_stiffness[0,1] * wgp * detJ
                k_e[I*2+1, J*2+0] += nodal_stiffness[1,0] * wgp * detJ
                k_e[I*2+1, J*2+1] += nodal_stiffness[1,1] * wgp * detJ


    return r_e, k_e


def Elmt_Post(XL, UL, Hn, Ht, Mat, dt, PostName):
    '''
    '''
    ## NEW post strategy is to return indexes
    # and contribuzions
    # 4 nodes
    r_post = np.zeros(4)

    # geometry and dofs
    XI = XL.reshape(-1,2)
    uI = UL.reshape(-1,2)
    

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
        r_post = np.array([UL[0], UL[2], UL[4], UL[6]])
        return r_post
    elif (PostName=="UY"):
        r_post = np.array([UL[1], UL[3], UL[5], UL[7]])
        return r_post
    elif (PostName=="SigMises"):
        return r_post
    else:
        print("Waring: PostName "+PostName+" not defined!")
        return np.array([0.0, 0.0, 0.0, 0.0])
