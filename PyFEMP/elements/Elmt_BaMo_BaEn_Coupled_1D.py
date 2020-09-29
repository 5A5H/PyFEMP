#from numpy import *
import numpy as np


def Elmt_Init():
    NoElementDim = 1
    NoElementNodes = 2
    NoElementHistory = 12
    ElementDofNames = ["U", "T"]
    ElementMaterialNames = ["E", "A", "nu", "rho", "alpha", "c"]
    ElementPostNames = ["A", "Sig", "q"]
    return NoElementDim, NoElementNodes, ElementDofNames, NoElementHistory, ElementMaterialNames, ElementPostNames


def Newmark_V(U, U_n, V_n, A_n, gamma, beta, dt):
    return gamma/(beta*dt) * (U-U_n) + (1.0-gamma/beta) * V_n + dt * (1.0-gamma/(2.0*beta)) * A_n


def Newmark_A(U, U_n, V_n, A_n, gamma, beta, dt):
    return 1.0/(beta*dt*dt) * (U-U_n - dt*V_n - dt*dt*(1.0/2.0 - beta)*A_n)


def Elmt_KS(XI, UI, Hn, Ht, Mat, dt):
    '''
    '''
    verbose = False  # True;
    # element vector /matrix
    r_e = np.zeros(4)
    k_e = np.zeros((4, 4))
    # nodal coordinates
    X1 = XI[0]
    X2 = XI[1]
    if (verbose):
        print('---------------')
        print('X1: ', X1)
        print('X2: ', X2)
    # nodal degrees of freedom
    U1 = UI[0]
    U2 = UI[2]
    T1 = UI[1]
    T2 = UI[3]
    # load previous values from history Hn
    U1_n = Hn[0]
    U2_n = Hn[1]
    T1_n = Hn[2]
    T2_n = Hn[3]
    V1_n = Hn[4]
    V2_n = Hn[5]
    T1dt_n = Hn[6]
    T2dt_n = Hn[7]
    A1_n = Hn[8]
    A2_n = Hn[9]
    T1ddt_n = Hn[10]
    T2ddt_n = Hn[11]
    if (verbose):
        print('U_n :', [U1_n, U2_n])
    if (verbose):
        print('V_n :', [V1_n, V2_n])
    if (verbose):
        print('A_n :', [A1_n, A2_n])
    # material parameter
    Emod = Mat[0]
    Area = Mat[1]
    mue = Mat[2]
    rho = Mat[3]
    alpha = Mat[4]
    c = Mat[5]
    # Newmark Time integration
    gamma = 1.0/2.0
    beta = 1.0/4.0
    V1 = Newmark_V(U1, U1_n, V1_n, A1_n, gamma, beta, dt)
    V2 = Newmark_V(U2, U2_n, V2_n, A2_n, gamma, beta, dt)
    T1dt = Newmark_V(T1, T1_n, T1dt_n, T1ddt_n, gamma, beta, dt)
    T2dt = Newmark_V(T2, T2_n, T1dt_n, T1ddt_n, gamma, beta, dt)
    A1 = Newmark_A(U1, U1_n, V1_n, A1_n, gamma, beta, dt)
    A2 = Newmark_A(U2, U2_n, V2_n, A2_n, gamma, beta, dt)
    T1ddt = Newmark_A(T1, T1_n, T1dt_n, T1ddt_n, gamma, beta, dt)
    T2ddt = Newmark_A(T2, T2_n, T2dt_n, T2ddt_n, gamma, beta, dt)
    if (verbose):
        print('U :', [U1, U2])
    if (verbose):
        print('V :', [V1, V2])
    if (verbose):
        print('A :', [A1, A2])
    # export for next time step
    Ht[0] = U1
    Ht[1] = U2
    Ht[2] = T1
    Ht[3] = T2
    Ht[4] = V1
    Ht[5] = V2
    Ht[6] = T1dt
    Ht[7] = T2dt
    Ht[8] = A1
    Ht[9] = A2
    Ht[10] = T1ddt
    Ht[11] = T2ddt
    # Geometry and mapping
    L = X2 - X1
    dxidx = 2.0/L
    # Shape functions, evaluated at xi = 0 (midpoint integration)
    N = 1.0/2.0 * np.array([1.0, 1.0])
    B = 1.0/2.0 * np.array([-1.0, 1.0]) * dxidx
    # Element quantities
    A = N[0] * A1 + N[1] * A2
    Eps = B[0] * U1 + B[1] * U2
    Eps_dt = B[0] * V1 + B[1] * V2
    T_dt = N[0] * T1dt + N[1] * T2dt
    Grad_T = B[0] * T1 + B[1] * T2
    Sig = Emod * Eps + mue * Eps_dt
    # Element vector
    r_e[0] += L * Area * (Sig * B[0] + rho * A * N[0])
    r_e[2] += L * Area * (Sig * B[1] + rho * A * N[1])
    # Element matrix
    k_e[0][0] += L * Area * Emod*B[0]*B[0] + L * Area * mue * \
        (gamma/(beta*dt))*B[0]*B[0] + L * \
        Area * rho*(1.0/(beta*dt**2))*N[0]*N[0]
    k_e[0][2] += L * Area * Emod*B[1]*B[0] + L * Area * mue * \
        (gamma/(beta*dt))*B[1]*B[0] + L * \
        Area * rho*(1.0/(beta*dt**2))*N[1]*N[0]
    k_e[2][0] += L * Area * Emod*B[0]*B[1] + L * Area * mue * \
        (gamma/(beta*dt))*B[0]*B[1] + L * \
        Area * rho*(1.0/(beta*dt**2))*N[0]*N[1]
    k_e[2][2] += L * Area * Emod*B[1]*B[1] + L * Area * mue * \
        (gamma/(beta*dt))*B[1]*B[1] + L * \
        Area * rho*(1.0/(beta*dt**2))*N[1]*N[1]
    # BaEn static term
    r_e[1] += L * Area * alpha * Grad_T * B[0]
    r_e[3] += L * Area * alpha * Grad_T * B[1]

    k_e[1][1] += L * Area * alpha * B[0] * B[0]
    k_e[1][3] += L * Area * alpha * B[0] * B[1]
    k_e[3][1] += L * Area * alpha * B[1] * B[0]
    k_e[3][3] += L * Area * alpha * B[1] * B[1]

    # BaEn dynamic term
    r_e[1] += L * Area * c * T_dt * N[0]
    r_e[3] += L * Area * c * T_dt * N[1]

    k_e[1][1] += L * Area * c * N[0] * N[0] * (gamma/(beta*dt))
    k_e[1][3] += L * Area * c * N[0] * N[1] * (gamma/(beta*dt))
    k_e[3][1] += L * Area * c * N[1] * N[0] * (gamma/(beta*dt))
    k_e[3][3] += L * Area * c * N[1] * N[1] * (gamma/(beta*dt))

    if (verbose):
        print('r_e : ', r_e)
    if (verbose):
        print('k_e : ')
    if (verbose):
        print(k_e)
    return r_e, k_e


def Elmt_Post(XI, UI, Hn, Ht, Mat, dt, PostName):
    '''
    '''
    # nodal coordinates
    X1 = XI[0]
    X2 = XI[1]
    # nodal degrees of freedom
    U1 = UI[0]
    U2 = UI[2]
    T1 = UI[1]
    T2 = UI[3]
    # load previous values from history Hn
    U1_n = Hn[0]
    U2_n = Hn[1]
    T1_n = Hn[2]
    T2_n = Hn[3]
    V1_n = Hn[4]
    V2_n = Hn[5]
    T1dt_n = Hn[6]
    T2dt_n = Hn[7]
    A1_n = Hn[8]
    A2_n = Hn[9]
    T1ddt_n = Hn[10]
    T2ddt_n = Hn[11]
    # material parameter
    Emod = Mat[0]
    Area = Mat[1]
    mue = Mat[2]
    rho = Mat[3]
    alpha = Mat[4]
    c = Mat[5]
    # Newmark Time integration
    gamma = 1.0/2.0
    beta = 1.0/4.0
    V1 = Newmark_V(U1, U1_n, V1_n, A1_n, gamma, beta, dt)
    V2 = Newmark_V(U2, U2_n, V2_n, A2_n, gamma, beta, dt)
    T1dt = Newmark_V(T1, T1_n, T1dt_n, T1ddt_n, gamma, beta, dt)
    T2dt = Newmark_V(T2, T2_n, T1dt_n, T1ddt_n, gamma, beta, dt)
    A1 = Newmark_A(U1, U1_n, V1_n, A1_n, gamma, beta, dt)
    A2 = Newmark_A(U2, U2_n, V2_n, A2_n, gamma, beta, dt)
    T1ddt = Newmark_A(T1, T1_n, T1dt_n, T1ddt_n, gamma, beta, dt)
    T2ddt = Newmark_A(T2, T2_n, T2dt_n, T2ddt_n, gamma, beta, dt)
    # export for next time step
    Ht[0] = U1
    Ht[1] = U2
    Ht[2] = T1
    Ht[3] = T2
    Ht[4] = V1
    Ht[5] = V2
    Ht[6] = T1dt
    Ht[7] = T2dt
    Ht[8] = A1
    Ht[9] = A2
    Ht[10] = T1ddt
    Ht[11] = T2ddt
    # Geometry and mapping
    L = X2 - X1
    dxidx = 2.0/L
    # Shape functions, evaluated at xi = 0 (midpoint integration)
    N = 1.0/2.0 * np.array([1.0, 1.0])
    B = 1.0/2.0 * np.array([-1.0, 1.0]) * dxidx
    # Element quantities
    A = N[0] * A1 + N[1] * A2
    Eps = B[0] * U1 + B[1] * U2
    Eps_dt = B[0] * V1 + B[1] * V2
    T_dt = N[0] * T1dt + N[1] * T2dt
    Sig = Emod * Eps + mue * Eps_dt
    Grad_T = B[0] * T1 + B[1] * T2
    q = - alpha * Grad_T
    if PostName == "Sig":
        return np.array([Sig, Sig])
    if PostName == "A":
        return np.array([A1, A2])
    if PostName == "q":
        return np.array([q, q])

    return np.array([0.0, 0.0])
