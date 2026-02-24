import numpy as np
from numba import njit
from TRILOC6 import TRILOC6   # must be njit compiled


@njit(cache=True)
def LOCIN6(XJ1, XJ2, XJ3, NODO, NCONEC, hd):

    HW = np.zeros(NCONEC, dtype=np.complex128)
    GW = np.zeros(NCONEC, dtype=np.complex128)
    HS = np.zeros(NCONEC, dtype=np.complex128)

    GII1 = np.zeros(15)
    GII2 = np.zeros(15)

    HDIF = 0.0

    # ----------------------------------
    # Natural coordinates (9-node quad)

    GI1 = np.array([-1.0, 0.0, 1.0, 1.0, 1.0, 0.0, -1.0, -1.0, 0.0])
    GI2 = np.array([-1.0,-1.0,-1.0, 0.0, 1.0, 1.0,  1.0,  0.0, 0.0])

    for i in range(8):
        GII1[i] = GI1[i]
        GII2[i] = GI2[i]

    for i in range(7):
        GII1[8+i] = GI1[i]
        GII2[8+i] = GI2[i]

    # ----------------------------------
    # Collocation node

    V1 = np.zeros(3)
    V2 = np.zeros(3)

    V1[0] = GI1[NODO]
    V2[0] = GI2[NODO]

    XI = np.zeros(3)
    XI[0] = XJ1[NODO]
    XI[1] = XJ2[NODO]
    XI[2] = XJ3[NODO]

    # ----------------------------------
    # Corner nodes (0,2,4,6)

    if (NODO == 0) or (NODO == 2) or (NODO == 4) or (NODO == 6):

        AREA = 2.0

        for k in range(1,3):
            M = 2*k

            V1[1] = GII1[NODO + M]
            V2[1] = GII2[NODO + M]

            V1[2] = GII1[NODO + M + 2]
            V2[2] = GII2[NODO + M + 2]

            HW, GW, HS = TRILOC6(
                HW, GW, HS,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2, NCONEC, hd
            )

    # ----------------------------------
    # Mid-side nodes (1,3,5,7)

    elif (NODO == 1) or (NODO == 3) or (NODO == 5) or (NODO == 7):

        for k in range(1,4):

            M = 2*k - 1

            V1[1] = GII1[NODO + M]
            V2[1] = GII2[NODO + M]

            V1[2] = GII1[NODO + M + 2]
            V2[2] = GII2[NODO + M + 2]

            AREA = 0.5 * (
                V1[1]*V2[2] + V1[0]*V2[1] + V1[2]*V2[0]
              - V1[1]*V2[0] - V1[2]*V2[1] - V1[0]*V2[2]
            )

            HW, GW, HS = TRILOC6(
                HW, GW, HS,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2, NCONEC, hd
            )

    # ----------------------------------
    # Center node (8)

    elif NODO == 8:

        AREA = 1.0

        for k in range(1,5):

            M = 2*k - 1

            V1[1] = GII1[M]
            V2[1] = GII2[M]

            V1[2] = GII1[M + 2]
            V2[2] = GII2[M + 2]

            HW, GW, HS = TRILOC6(
                HW, GW, HS,
                XI, XJ1, XJ2, XJ3,
                AREA, V1, V2, NCONEC, hd
            )

    # ----------------------------------

    HDIF = HDIF + HW[NODO] - HS[NODO]

    return HW, GW, HS, HDIF
