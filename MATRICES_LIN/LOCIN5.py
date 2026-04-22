import numpy as np
from numba import njit

from .TRILOC5 import TRILOC5

@njit(fastmath=True)
def LOCIN5(XJ1, XJ2, XJ3, NODO, NCONEC, DEPTH):

    # Element integrals
    HW = np.zeros(NCONEC, dtype=np.complex128)
    GW = np.zeros(NCONEC, dtype=np.complex128)

    # Homogeneous coordinates vectors
    GI1 = np.zeros(NCONEC)
    GI2 = np.zeros(NCONEC)
    
    GII1 = np.zeros(15)
    GII2 = np.zeros(15)
    
    # Homogeneous coordinates
    GI1[0] = -1.0
    GI1[1] =  1.0
    GI1[2] =  1.0
    GI1[3] = -1.0

    GI2[0] = -1.0
    GI2[1] = -1.0
    GI2[2] =  1.0
    GI2[3] =  1.0

    # Copy nodes 1..8
    for i in range(4):
        GII1[i] = GI1[i]
        GII2[i] = GI2[i]

    # Copy nodes 1..7 again into 9..15
    for i in range(3):
        GII1[4+i] = GI1[i]
        GII2[4+i] = GI2[i]

    # Triangle vertices
    V1 = np.zeros(3)
    V2 = np.zeros(3)

    # First triangle vertex
    V1[0] = GI1[NODO]
    V2[0] = GI2[NODO]

    # Field point
    XI = np.zeros(3)
    
    XI[0] = XJ1[NODO]
    XI[1] = XJ2[NODO]
    XI[2] = XJ3[NODO]

    # ---- Subdivision cases ----
    # Corner nodes

    AREA = 2.0

    for k in range(2):
        
        M = k + 1

        V1[1] = GII1[NODO + M]
        V2[1] = GII2[NODO + M]

        V1[2] = GII1[NODO + M + 1]
        V2[2] = GII2[NODO + M + 1]

        HW, GW = TRILOC5(
            HW, GW,
            XI, XJ1, XJ2, XJ3,
            AREA, V1, V2,
            NCONEC, DEPTH
        )

    return HW, GW