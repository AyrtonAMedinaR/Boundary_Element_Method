import numpy as np
from numba import njit

@njit(fastmath=False)
def FUNDA6(XI, ETA, XG1, XG2, XG3, K, hd):

    # --- Field point ---
    XP = XI[0]
    YP = XI[1]
    ZP = XI[2]

    # --- Source (Gauss) point ---
    XG = XG1[K]
    YG = XG2[K]
    ZG = XG3[K]

    # Small tolerance
    eps = 1e-12

    # ---------------------------
    # Primary image
    RD1 = XG - XP
    RD2 = YG - YP
    RD3 = ZG - ZP

    R2 = RD1*RD1 + RD2*RD2 + RD3*RD3
    # if R2 < eps:
    #     R2 = eps
    R = np.sqrt(R2)

    RDN = (RD1*ETA[0] + RD2*ETA[1] + RD3*ETA[2]) / R

    # ---------------------------
    # Mirror image
    RD1_2 = XG - XP
    RD2_2 = YG - YP
    RD3_2 = ZG + ZP + 2.0*abs(hd)

    R2_2 = RD1_2*RD1_2 + RD2_2*RD2_2 + RD3_2*RD3_2
    if R2_2 < eps:
        R2_2 = eps
    R_2 = np.sqrt(R2_2)

    RDN_2 = (RD1_2*ETA[0] + RD2_2*ETA[1] + RD3_2*ETA[2]) / R_2

    # ---------------------------
    # Constants
    c = 1.0 / (4.0*np.pi)

    # ---------------------------
    # Fundamental solutions
    U_real = c*(1.0/R + 1.0/R_2)

    Q_real = -c*((RDN/R)/R + (RDN_2/R_2)/R_2)

    QS = -c*(RDN/(R*R) + RDN_2/(R_2*R_2))

    # Force complex type
    U = U_real + 0j
    Q = Q_real + 0j

    return U, Q, QS



