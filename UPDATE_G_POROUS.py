import numpy as np

def UPDATE_G_POROUS(Coef, k, omega, eps, l, Gporous, NCONEC, R0_STRUCT, R1_STRUCT, POS_M, POS_N, N, Region_M, Region_N):
    
    N_prev_M  = sum(N[0:Region_M])
    N_prev_N  = sum(N[0:Region_N])

    varphi = 0.6 + 0.4 * eps**2
    alpha = (1/(eps*varphi) - 1)**2

    R_M = R0_STRUCT[:,0].reshape(-1, NCONEC)
    R_N = R1_STRUCT[:,0].reshape(-1, NCONEC)

    dPHI = 1j*k*Gporous*(Coef[R_M+N_prev_M] - Coef[R_N+N_prev_N])

    Sup = np.abs(dPHI)**3
    Inf = np.abs(dPHI)**2

    n1 = R_M[:,0]
    n2 = R_M[:,1]
    n8 = R_M[:,7]

    dx = np.linalg.norm(POS_M[n1] - POS_M[n2], axis=1)
    dy = np.linalg.norm(POS_M[n1] - POS_M[n8], axis=1)

    w = dx * dy / 4

    NOM = w * (
        Sup[:,0] + Sup[:,2] + Sup[:,4] + Sup[:,6]
        + 2*(Sup[:,1] + Sup[:,3] + Sup[:,5] + Sup[:,7])
        + 4*Sup[:,8]
    )

    DEN = w * (
        Inf[:,0] + Inf[:,2] + Inf[:,4] + Inf[:,6]
        + 2*(Inf[:,1] + Inf[:,3] + Inf[:,5] + Inf[:,7])
        + 4*Inf[:,8]
    )

    betta = (alpha/2) * (8/(3*np.pi)) * NOM.sum() / DEN.sum()
    R = betta * k / omega

    Gporous_new = 1 / (R - 1j*k*l)

    return 0.5 * (Gporous_new + Gporous)

