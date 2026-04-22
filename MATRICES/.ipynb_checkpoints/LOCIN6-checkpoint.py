import numpy as np
from numba import njit

from .TRILOC6 import TRILOC6

@njit(fastmath=True)
def LOCIN6(XJ1, XJ2, XJ3, NODO, NCONEC, DEPTH):

    # Element integrals
    HW = np.zeros(NCONEC, dtype=np.complex128)
    GW = np.zeros(NCONEC, dtype=np.complex128)

    # Homogeneous coordinates vectors
    GI1 = np.zeros(9)
    GI2 = np.zeros(9)
    
    GII1 = np.zeros(15)
    GII2 = np.zeros(15)
    
    # Homogeneous coordinates
    GI1[0] = -1.0; GI1[1] = 0.0; GI1[2] = 1.0
    GI1[3] = 1.0;  GI1[4] = 1.0; GI1[5] = 0.0
    GI1[6] = -1.0; GI1[7] = -1.0; GI1[8] = 0.0

    GI2[0] = -1.0; GI2[1] = -1.0; GI2[2] = -1.0
    GI2[3] = 0.0;  GI2[4] = 1.0;  GI2[5] = 1.0
    GI2[6] = 1.0;  GI2[7] = 0.0;  GI2[8] = 0.0

    # Copy nodes 1..8
    for i in range(8):
        GII1[i] = GI1[i]
        GII2[i] = GI2[i]

    # Copy nodes 1..7 again into 9..15
    for i in range(7):
        GII1[8+i] = GI1[i]
        GII2[8+i] = GI2[i]

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
    if NODO == 0 or NODO == 2 or NODO == 4 or NODO == 6:

        AREA = 2.0

        for k in range(2):

            M = 2*(k+1)

            V1[1] = GII1[NODO + M]
            V2[1] = GII2[NODO + M]

            V1[2] = GII1[NODO + M + 2]
            V2[2] = GII2[NODO + M + 2]

            HW, GW = TRILOC6(
                HW, GW,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2,
                NCONEC, DEPTH
            )

    # Edge nodes
    elif NODO == 1 or NODO == 3 or NODO == 5 or NODO == 7:

        for k in range(3):

            M = 2*(k+1) - 1

            V1[1] = GII1[NODO + M]
            V2[1] = GII2[NODO + M]

            V1[2] = GII1[NODO + M + 2]
            V2[2] = GII2[NODO + M + 2]

            AREA = 0.5 * (
                V1[1]*V2[2] + V1[0]*V2[1] + V1[2]*V2[0]
                - V1[1]*V2[0] - V1[2]*V2[1] - V1[0]*V2[2]
            )

            HW, GW = TRILOC6(
                HW, GW,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2,
                NCONEC, DEPTH
            )

    # Center node
    elif NODO == 8:

        AREA = 1.0

        for k in range(4):

            M = 2*(k+1) - 1
           
            V1[1] = GII1[M-1]
            V2[1] = GII2[M-1]

            V1[2] = GII1[M+1]
            V2[2] = GII2[M+1]

            HW, GW = TRILOC6(
                HW, GW,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2,
                NCONEC, DEPTH
            )

    return HW, GW