import numpy as np
from numba import njit

@njit(fastmath=True)
def FUNDA5(XP, YP, ZP, ETA, XG, YG, ZG, DEPTH):
    
    # Fundamental solution evaluated at a Gauss point (3D Laplace Green function)

    # Parameters
    # XP, YP, ZP : Source point coordinates (collocation point)
    # ETA : array-like (3,) Outward normal vector at the element
    # XG, YG, ZG : Gauss point coordinates

    # Returns
    # U : complex, Green function 1/(4*pi*R)
    # Q : complex, Normal derivative of Green function

    # Distance components
    RD1 = XG - XP
    RD2 = YG - YP
    RD3 = ZG - ZP

    # Distance
    R = np.sqrt(RD1**2 + RD2**2 + RD3**2)

    # Dot product divided by distance
    RDN = (RD1*ETA[0] + RD2*ETA[1] + RD3*ETA[2]) / R

    # ---- Image source ----
    RD1_2 = RD1
    RD2_2 = RD2
    RD3_2 = ZG + ZP + 2.0*DEPTH

    R_2 = np.sqrt(RD1_2**2 + RD2_2**2 + RD3_2**2)

    RDN_2 = (RD1_2*ETA[0] + RD2_2*ETA[1] + RD3_2*ETA[2]) / R_2

    # ---- Constants ----
    c = 1.0 / (4.0 * np.pi)

    invR = 1.0 / R
    invR2 = 1.0 / R_2

    # ---- Fundamental solutions ----
    U_real = c * (invR + invR2)

    Q_real = -c * (RDN * invR * invR + RDN_2 * invR2 * invR2)    

    # complex outputs (required by BEM matrix assembly)
    U = complex(U_real, 0.0)
    Q = complex(Q_real, 0.0)
    
    return U, Q