# T1 Plasti Element
import numpy as np


def Elmt_Init():
    NoElementDim = 2
    NoElementNodes = 3
    NoElementHistory = 7 
    ElementDofNames = ["UX", "UY"]
    ElementMaterialNames = ["E", "nu", "y0", "kh"]
    ElementPostNames = ["UX", "UY", "SigMises", "a"]
    return NoElementDim, NoElementNodes, ElementDofNames, NoElementHistory, ElementMaterialNames, ElementPostNames


def SH0_T1(xi, eta):
    '''
    SH0_T1(xi, eta) -> SH0
    Return a two dimensional array containing shape functions and derivatives for T1 element.
    Usage: SH0( NodeNumber, SHPIndex)  
    with SHPIndex = {
        0 -> shape function, 
        1 -> derived shape function w.r.t. xi, 
        2 -> derived shape function w.r.t. eta
        } 
    '''
    return np.array([
            [xi        ,  1.0,   0.0],
            [eta       ,  0.0,   1.0],
            [1 -xi -eta, -1.0,  -1.0],
        ], dtype=np.float64)

def BmatVoigt2D(SHP):
    '''
    BmatVoigt(SHP) -> Bmat
    Returns a B-Matrix (as dim:3) for computing the strain vector in voigt notation.
    This B-Matrix assumes a 2D plane strain approximation.
    SHP is a shape function matrix with derived functions w.r.t. physical space.
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

def HookeMatVoigt(kappa, mue):
    '''
    HookeMatVoigt(kappa, mue) -> Cmat
    Returns the constitutive Voigt MATRIX(6,6) for a Hooke material law.
    The input are the kompression modulus kappa and the shear modulus mue.
    sig_i = Cmat_ij * eps_j
    with
    sig_i = [sig_11, sig_22, sig_33, sig_12, sig_23, sig_13]
    eps_i = [eps_11, eps_22, eps_33, 2*eps_12, 2*eps_23, 2*eps_13]
    '''
    return np.array([
        [kappa+(4.0/3.0)*mue, kappa-(2.0/3.0)*mue, kappa-(2.0/3.0)*mue, 0  , 0  , 0  ],
        [kappa-(2.0/3.0)*mue, kappa+(4.0/3.0)*mue, kappa-(2.0/3.0)*mue, 0  , 0  , 0  ],
        [kappa-(2.0/3.0)*mue, kappa-(2.0/3.0)*mue, kappa+(4.0/3.0)*mue, 0  , 0  , 0  ],
        [0                  , 0                  , 0                  , mue, 0  , 0  ],
        [0                  , 0                  , 0                  , 0  , mue, 0  ],
        [0                  , 0                  , 0                  , 0  , 0  , mue]
    ], dtype=np.float64)


# definition of auxillary matrix
pp1, pp2, pp3 = (2.0/3.0), -(1.0/3.0), (1.0/2.0)
PP = np.array( \
[[pp1, pp2, pp2, 0.0, 0.0, 0.0], \
[ pp2, pp1, pp2, 0.0, 0.0, 0.0], \
[ pp2, pp2, pp1, 0.0, 0.0, 0.0], \
[ 0.0, 0.0, 0.0, pp3, 0.0, 0.0], \
[ 0.0, 0.0, 0.0, 0.0, pp3, 0.0], \
[ 0.0, 0.0, 0.0, 0.0, 0.0, pp3]] \
,dtype=np.float64)

def Elmt_KS(XL, UL, Hn, Ht, Mat, dt):
    '''
    '''
    verbose = False  # True;
    if verbose: print('XI :',XL)    # XL  = [x11, x12, x21, x22, x31, x32]
    if verbose: print('UI :',UL)    # UL  = [u11, u12, u21, u22, u31, u32]
    if verbose: print('Hn :',Hn)    # Hn  = [eps_pl_11, eps_pl_22, eps_pl_33, 2*eps_pl_12, 2*eps_pl_23, 2*eps_pl_13, a]
    if verbose: print('Ht :',Ht)    # Ht  = [eps_pl_11, eps_pl_22, eps_pl_33, 2*eps_pl_12, 2*eps_pl_23, 2*eps_pl_13, a]
    if verbose: print('b  :',Mat)   # Mat = [Emod, nue, y0, kh]
    if verbose: print('dt :',dt)    # dt  = dt
    
    # element specific paremeters
    NoElementNodes = 3
    NoNodalDOF     = 2
    NoDimension    = 2
    
    # initialize element vector /matrix
    r_e = np.zeros(NoElementNodes*NoNodalDOF)
    k_e = np.zeros((NoElementNodes*NoNodalDOF, NoElementNodes*NoNodalDOF))

    # geometry and dofs
    XI = XL.reshape(-1, NoDimension)
    uI = UL.reshape(-1, NoNodalDOF)

    # Material Parameters
    Emod, nu, y0, kh = Mat[0], Mat[1], Mat[2], Mat[3]
    kappa , mue = Emod/(3*(1.0-2.0*nu)), Emod/(2.0*(1.0+nu))

    # constitutive matrix (hooke)
    Cmat = HookeMatVoigt(kappa, mue)

    # provide integration points
    EGP = np.array([[(1.0/3.0), (1.0/3.0), (1.0/2.0)]])
    NoInt = len(EGP)

    # start integration Loop
    for GP in range(NoInt):
        if verbose: print('GP: ',GP)
        xi, eta, wgp  = EGP[GP]

        # read gp history
        NoHGP = 7 # number of history at each gp
        eps_pl_n = Hn[ GP*NoHGP : GP*NoHGP+6]
        a_n = Hn[ GP*NoHGP+6]

        # compute current shape functions
        SH0 = SH0_T1(xi, eta)

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
        if verbose: print(np.array([[eps[0], eps[3]/2, eps[4]/2], [eps[3]/2, eps[1], eps[5]/2], [eps[4]/2, eps[5]/2, eps[2]]]))

        ###############################################
        ############# begin plastic part ##############
        ###############################################
        
        # compute elastic trail stresses
        eps_el_tr = eps - eps_pl_n
        if verbose: print('eps_el_tr: ')
        if verbose: print(np.array([[eps_el_tr[0], eps_el_tr[3]/2, eps_el_tr[4]/2], [eps_el_tr[3]/2, eps_el_tr[1], eps_el_tr[5]/2], [eps_el_tr[4]/2, eps_el_tr[5]/2, eps_el_tr[2]]]))

        # compute deviatoric trail stresses
        sig_tr = np.einsum('ij,j->i'  , Cmat, eps_el_tr)
        tr_sig_tr  = sig_tr[0] + sig_tr[1] + sig_tr[2]
        dev_sig_tr = sig_tr - (1.0/3.0) * tr_sig_tr * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

        # compute norm of deviatoric trail stresses
        norm_dev_sig_tr = np.sqrt(
                  dev_sig_tr[0]**2 +       dev_sig_tr[1]**2 +       dev_sig_tr[2]**2 + \
            2.0 * dev_sig_tr[3]**2 + 2.0 * dev_sig_tr[4]**2 + 2.0 * dev_sig_tr[5]**2
        )

        # compute yield criterion
        phi_tr = norm_dev_sig_tr - np.sqrt(2.0/3.0) * (y0 + (2.0/3.0)*kh*a_n)
        if verbose: print('norm_dev_sig_tr: ', norm_dev_sig_tr)
        if verbose: print('yield:           ', np.sqrt(2.0/3.0) * (y0 + (2.0/3.0)*kh*a_n))
        if verbose: print('phi_tr:          ', phi_tr)

        # check yield criterion
        if (phi_tr > 1e-8): # elasto-plastic loading
            
            # compute plastic strain increment
            delta_a = phi_tr/(2.0*mue + (2.0/3.0*kh))
            a = a_n + np.sqrt(2.0/3.0) * delta_a
            if verbose: print('delta_a:         ', delta_a)

            # compute plastic flow director
            n_tr = dev_sig_tr/norm_dev_sig_tr
            if verbose: print('n_tr: ')
            if verbose: print(np.array([[n_tr[0], n_tr[3], n_tr[4]], [n_tr[3], n_tr[1], n_tr[5]], [n_tr[4], n_tr[5], n_tr[2]]]))

            # compute plastic strain - take care of voigt notation!
            eps_pl = eps_pl_n + delta_a * n_tr * np.array([1.0, 1.0, 1.0, 2.0, 2.0, 2.0]) 
            if verbose: print('eps_pl: ')
            if verbose: print(np.array([[eps_pl[0], eps_pl[3]/2, eps_pl[4]/2], [eps_pl[3]/2, eps_pl[1], eps_pl[5]/2], [eps_pl[4]/2, eps_pl[5]/2, eps_pl[2]]]))

            # compute new elastic strains
            eps_el = eps - eps_pl
            if verbose: print('eps_el: ')
            if verbose: print(np.array([[eps_el[0], eps_el[3]/2, eps_el[4]/2], [eps_el[3]/2, eps_el[1], eps_el[5]/2], [eps_el[4]/2, eps_el[5]/2, eps_el[2]]]))

            sig = np.einsum('ij,j->i', Cmat, eps_el)

            # modification of material tangent
            fact1 = 1.0 - (2.0*mue*delta_a/norm_dev_sig_tr)
            fact2 = 2.0*mue/(2.0*mue + (2.0/3.0*kh)) - (2.0*mue*delta_a/norm_dev_sig_tr)
            Cmat += -2.0*mue*PP \
                    + 2.0*mue*PP*fact1 \
                    - 2.0*mue*np.einsum('i,j->ij', n_tr, n_tr)*fact2

            # export history
            Ht[ GP*NoHGP : GP*NoHGP+6] = eps_pl
            Ht[ GP*NoHGP+6] = a

        else: # elastic loading

            # old plastic state is stil valid
            eps_pl = eps_pl_n
            a = a_n

            # compute elastic strains
            eps_el = eps-eps_pl

            sig = np.einsum('ij,j->i', Cmat, eps_el)

        ###############################################
        ############## end plastic part ###############
        ############################################### 

        # export right hand side | this element has 3 nodes with 2 dofs each
        for I in range(NoElementNodes):
            
            # compute nodal right hand side
            nodal_rhs_vec    = np.einsum('i,ij->j',sig, Bmat[I])

            # integrate nodal right hand side and export
            r_e[I*2+0] += nodal_rhs_vec[0] * wgp * detJ
            r_e[I*2+1] += nodal_rhs_vec[1] * wgp * detJ

            for J in range(NoElementNodes):

                # compute nodal stiffness matrix
                nodal_stiffness = np.einsum('ki,ko,oj->ij', Bmat[I], Cmat, Bmat[J])

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
    # 3 nodes
    r_post = np.zeros(3)

    # geometry and dofs
    XI = np.array([[XL[0], XL[1]], [XL[2], XL[3]], [XL[ 4], XL[ 5]]], dtype=np.float64)
    uI = np.array([[UL[0], UL[1]], [UL[2], UL[3]], [UL[ 4], UL[ 5]]], dtype=np.float64)

    # Material Parameters
    Emod, nu, y0, kh = Mat[0], Mat[1], Mat[2], Mat[3]
    lam , mue = (Emod*nu)/((1.0+nu)*(1.0-2.0*nu)), Emod/(2.0*(1.0+nu))
    kappa = lam + (2.0/3.0)*mue

    # constitutive matrix (hooke)
    Cmat = HookeMatVoigt(kappa, mue)

    # provide integration points
    EGP = np.array([[(1.0/3.0), (1.0/3.0), (1.0/2.0)]])
    NoInt = len(EGP)

    # start integration Loop
    for GP in range(NoInt):
        xi, eta, wgp  = EGP[GP]

        # read gp history
        NoHGP = 7 # number of history at each gp
        eps_pl = Ht[ GP*NoHGP : GP*NoHGP+6]
        a = Ht[ GP*NoHGP+6]

        # compute current shape functions
        SH0 = SH0_T1(xi, eta)

        # compute mapping
        Jed = np.einsum('Ii,Ij->ij', XI ,SH0[:,1:3])
        detJ = np.linalg.det(Jed)
        if (detJ <= 0): raise NameError("Error unphysical mapping detected.")
        Jed_inv = np.linalg.inv(Jed)
        
        # map shape function derivative
        SHP = np.copy(SH0)
        SHP[:,1:3] = np.einsum('Ij,ji->Ii', SH0[:,1:3], Jed_inv)
        Bmat = BmatVoigt2D(SHP)

        # compute strains
        eps = np.einsum('Iij,Ij->i', Bmat, uI)

        # compute elastic strain
        eps_el = eps - eps_pl

        # compute stresses
        sig = np.einsum('ij,j->i'  , Cmat, eps_el)
        tr_sig  = sig[0] + sig[1] + sig[2]
        dev_sig = sig - (1.0/3.0) * tr_sig * np.array([1.0, 1.0, 1.0, 0.0, 0.0, 0.0])

        sig_vm = np.sqrt( \
            sig[0]**2 + sig[1]**2 + sig[2]**2 \
            -sig[0]*sig[1] -sig[0]*sig[2] -sig[1]*sig[2] \
            + 3* (sig[3]**2 + sig[4]**2 + sig[5]**2) \
            )
        
        if (PostName=="SigMises"):
            r_post += sig_vm

        if (PostName=="a"):
            r_post += a

    if (PostName=="UX"):
        r_post = np.array([UL[0], UL[2], UL[4]])
        return r_post
    elif (PostName=="UY"):
        r_post = np.array([UL[1], UL[3], UL[5]])
        return r_post
    elif (PostName=="SigMises"):
        return r_post
    elif (PostName=="a"):
        return r_post
    else:
        print("Waring: PostName "+PostName+" not defined!")
        return np.array([0.0, 0.0, 0.0])
