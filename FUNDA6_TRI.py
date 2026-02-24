import numpy as np
from numba import njit

@njit(fastmath=False)
def FUNDA6_TRI(XP, YP, ZP,
                  ET1, ET2, ET3,
                  XG, YG, ZG,
                  hd):

    eps = 1e-12

    # -------- Primary distance
    dx = XG - XP
    dy = YG - YP
    dz = ZG - ZP

    R2 = dx*dx + dy*dy + dz*dz
    if R2 < eps:
        R2 = eps
    R = np.sqrt(R2)

    invR = 1.0 / R

    RDN = (dx*ET1 + dy*ET2 + dz*ET3) * invR

    # -------- Image distance
    dz2 = ZG + ZP + 2.0*abs(hd)

    R2_2 = dx*dx + dy*dy + dz2*dz2
    if R2_2 < eps:
        R2_2 = eps
    R_2 = np.sqrt(R2_2)

    invR2 = 1.0 / R_2

    RDN_2 = (dx*ET1 + dy*ET2 + dz2*ET3) * invR2

    inv4pi = 1.0 / (4.0*np.pi)

    # -------- Fundamental solutions
    U_real = inv4pi * (invR + invR2)

    Q_real = -inv4pi * (RDN*invR*invR + RDN_2*invR2*invR2)

    QS = -inv4pi * (RDN/R2 + RDN_2/R2_2)

    # Force complex
    U = U_real + 0j
    Q = Q_real + 0j

    return U, Q, QS
