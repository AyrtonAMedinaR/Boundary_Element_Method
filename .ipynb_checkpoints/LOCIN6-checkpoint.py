import numpy as np
from TRILOC6 import TRILOC6   # must be njit compiled

def LOCIN6(XJ1, XJ2, XJ3, NODO, NCONEC, hd):

    # Convert MATLAB node index (1..9) to Python (0..8)
    nodo = NODO

    HW = np.zeros(NCONEC, dtype=complex)
    GW = np.zeros(NCONEC, dtype=complex)
    HS = np.zeros(NCONEC)

    GII1 = np.zeros(15)
    GII2 = np.zeros(15)

    HDIF = 0.0 + 0.0j

    GI1 = np.zeros(9)
    GI2 = np.zeros(9)

    # Homogeneous coordinates
    GI1[:] = [-1.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, 0.0]
    GI2[:] = [-1.0,-1.0,-1.0, 0.0, 1.0, 1.0,  1.0,  0.0, 0.0]

    # Copy nodes 1..8
    for i in range(8):
        GII1[i] = GI1[i]
        GII2[i] = GI2[i]

    # Copy nodes 1..7 again into 9..15
    for i in range(7):
        GII1[8+i] = GI1[i]
        GII2[8+i] = GI2[i]

    # Collocation node homogeneous coordinates
    V1 = np.zeros(3)
    V2 = np.zeros(3)
    XI = np.zeros(3)

    V1[0] = GI1[nodo]
    V2[0] = GI2[nodo]

    XI[0] = XJ1[nodo]
    XI[1] = XJ2[nodo]
    XI[2] = XJ3[nodo]

    # ---- Subdivision cases ----

    # Corner nodes
    if nodo in [0,2,4,6]:

        AREA = 2.0

        for k in range(2):
            M = 2*(k+1)

            V1[1] = GII1[nodo + M]
            V2[1] = GII2[nodo + M]

            V1[2] = GII1[nodo + M + 2]
            V2[2] = GII2[nodo + M + 2]

            HW, GW, HS = TRILOC6(
                HW, GW, HS,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2,
                NCONEC, hd
            )

    # Edge nodes
    elif nodo in [1,3,5,7]:

        for k in range(3):

            M = 2*(k+1) - 1

            V1[1] = GII1[nodo + M]
            V2[1] = GII2[nodo + M]

            V1[2] = GII1[nodo + M + 2]
            V2[2] = GII2[nodo + M + 2]

            AREA = 0.5 * (
                V1[1]*V2[2] + V1[0]*V2[1] + V1[2]*V2[0]
                - V1[1]*V2[0] - V1[2]*V2[1] - V1[0]*V2[2]
            )

            HW, GW, HS = TRILOC6(
                HW, GW, HS,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2,
                NCONEC, hd
            )

    # Center node
    elif nodo == 8:

        AREA = 1.0

        for k in range(4):

            M = 2*(k+1) - 1

            V1[1] = GII1[M-1]
            V2[1] = GII2[M-1]

            V1[2] = GII1[M+1]
            V2[2] = GII2[M+1]

            HW, GW, HS = TRILOC6(
                HW, GW, HS,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2,
                NCONEC, hd
            )

    # ---- Final correction ----
    HDIF = HDIF + HW[nodo] - HS[nodo]

    return HW, GW, HS, HDIF
