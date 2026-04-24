import numpy as np
from numba import njit
import math

@njit(fastmath=True)
def FUNDA5(XI, ETA, X, DEPTH):

    # --- Field (collocation) point ---
    XP = XI[0]
    YP = XI[1]
    ZP = XI[2]

    # --- Source (Gauss) point ---
    XG = X[0]
    YG = X[1]
    ZG = X[2]    

    # ---- Direct distance ----
    RD1 = XG - XP
    RD2 = YG - YP
    RD3 = ZG - ZP

    R = math.sqrt(RD1*RD1 + RD2*RD2 + RD3*RD3)

    RDN = (RD1*ETA[0] + RD2*ETA[1] + RD3*ETA[2]) / R

    # ---- Image source ----
    RD1_2 = RD1
    RD2_2 = RD2
    RD3_2 = ZG + ZP + 2.0*DEPTH

    R_2 = math.sqrt(RD1_2*RD1_2 + RD2_2*RD2_2 + RD3_2*RD3_2)

    RDN_2 = (RD1_2*ETA[0] + RD2_2*ETA[1] + RD3_2*ETA[2]) / R_2

    # ---- Constants ----
    c = 1.0 / (4.0 * math.pi)

    invR = 1.0 / R
    invR2 = 1.0 / R_2

    # ---- Fundamental solutions ----
    U_real = c * (invR + invR2)

    Q_real = -c * (RDN * invR * invR + RDN_2 * invR2 * invR2)    

    # complex outputs (required by BEM matrix assembly)
    U = complex(U_real, 0.0)
    Q = complex(Q_real, 0.0)

    return U, Q