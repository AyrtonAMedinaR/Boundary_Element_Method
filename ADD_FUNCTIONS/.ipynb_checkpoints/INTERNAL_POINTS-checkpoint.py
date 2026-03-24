import numpy as np

from MATRICES.ETAXJA import ETAXJA
from MATRICES.EXTIN6 import EXTIN6

def INTERNAL_POINTS(FI, DFI, hd, CX, CY, CZ, POS, NCONEC, KCONEC, NE_total, Region):

    NE = NE_total[Region]
    
    X = POS[:, 0]
    Y = POS[:, 1]
    Z = POS[:, 2]    

    L = len(CX)

    # Initialize DSOL
    DSOL = np.zeros(L, dtype=complex)

    # Loop over elements
    for IE in range(NE):

        # Allocate local arrays
        XJ1 = np.zeros(NCONEC)
        XJ2 = np.zeros(NCONEC)
        XJ3 = np.zeros(NCONEC)

        # Extract element node coordinates
        for II in range(NCONEC):
            IK = KCONEC[II, IE]  # already zero-based
            XJ1[II] = X[IK]
            XJ2[II] = Y[IK]
            XJ3[II] = Z[IK]

        # Call ETAXJA
        ETA1, ETA2, ETA3, XJA, XG1, XG2, XG3, FF = ETAXJA(
            XJ1, XJ2, XJ3, NCONEC
        )

        # Loop over internal points
        for k in range(L):
            XI = np.array([CX[k], CY[k], CZ[k]])

            # Call EXTIN6
            HW, GW, HS = EXTIN6(
                ETA1, ETA2, ETA3, XJA,
                XG1, XG2, XG3, FF,
                XI, NCONEC, hd
            )

            # Accumulate solution
            for M in range(NCONEC):
                JK = KCONEC[M, IE]
                DSOL[k] += GW[M] * DFI[M, IE] - HW[M] * FI[JK]

    return DSOL


