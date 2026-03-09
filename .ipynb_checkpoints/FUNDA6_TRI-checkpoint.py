import numpy as np

def FUNDA6_TRI(XI, ETA, X, hd):

    XP = XI[0]
    YP = XI[1]
    ZP = XI[2]

    XG = X[0]
    YG = X[1]
    ZG = X[2]

    # ---- Direct distance ----
    RD1 = XG - XP
    RD2 = YG - YP
    RD3 = ZG - ZP

    R = np.sqrt(RD1**2 + RD2**2 + RD3**2)

    RDN = (RD1*ETA[0] + RD2*ETA[1] + RD3*ETA[2]) / R

    # ---- Image source (method of images) ----
    RD1_2 = XG - XP
    RD2_2 = YG - YP
    RD3_2 = ZG + ZP + 2*np.abs(hd)

    R_2 = np.sqrt(RD1_2**2 + RD2_2**2 + RD3_2**2)

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
