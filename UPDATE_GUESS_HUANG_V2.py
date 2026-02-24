import numpy as np

def UPDATE_GUESS_HUANG_V2(Coef, N0, N01,
                               STORE_NODE_R1_POS, STORE_NODE_R2_POS,
                               k, omega, eps, l, Gporous,
                               POS1, NCONEC, POS2):

    varphi = 0.6 + 0.4 * eps**2
    alpha = (1/(eps*varphi) - 1)**2

    R1 = STORE_NODE_R1_POS[:,0].reshape(-1, NCONEC)
    R2 = STORE_NODE_R2_POS[:,0].reshape(-1, NCONEC)

    dPHI = 1j*k*Gporous*(Coef[R1+N0] - Coef[R2+N01])

    Sup = np.abs(dPHI)**3
    Inf = np.abs(dPHI)**2

    n1 = R1[:,0]
    n2 = R1[:,1]
    n8 = R1[:,7]

    dx = np.linalg.norm(POS1[n1] - POS1[n2], axis=1)
    dy = np.linalg.norm(POS1[n1] - POS1[n8], axis=1)

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

