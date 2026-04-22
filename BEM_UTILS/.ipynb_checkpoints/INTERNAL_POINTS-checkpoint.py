import numpy as np

from MATRICES.ETAXJA import ETAXJA
from MATRICES.EXTIN6 import EXTIN6

def INTERNAL_POINTS(FI, DFI, hd, CX, CY, CZ, POS, NCONEC, KCONEC, N, NE, Region):

    # TAKE THE PHI AND dPHI/dn THAT WILL BE USED IN A SPECIFIC REGION
    N_prev = sum(N[0:Region])
    N_act  = sum(N[0:Region + 1])
    FI_act = FI[N_prev:N_act]
    
    NE_prev = sum(NE[0:Region])
    NE_act  = sum(NE[0:Region + 1])
    DFI_act = DFI[:,NE_prev:NE_act]
    
    # TAKE THE KCONEC THAT WILL BE USED IN A SPECIFIC REGION
    KCONEC_act = KCONEC # [:,NE_prev:NE_act] - N_prev

    # TAKE THE NUMBER OF ELEMENTS NE THAT WILL BE USED IN A SPECIFIC REGION
    NE_act = NE[Region]
    
    X = POS[:, 0]
    Y = POS[:, 1]
    Z = POS[:, 2]    

    L = len(CX)

    # Initialize DSOL
    DSOL = np.zeros(L, dtype=complex)

    # Loop over elements
    for IE in range(NE_act):

        # Allocate local arrays
        XJ1 = np.zeros(NCONEC)
        XJ2 = np.zeros(NCONEC)
        XJ3 = np.zeros(NCONEC)

        # Extract element node coordinates
        for II in range(NCONEC):
            IK = KCONEC_act[II, IE]  # already zero-based
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
            HW, GW= EXTIN6(
                ETA1, ETA2, ETA3, XJA,
                XG1, XG2, XG3, FF,
                XI, NCONEC, hd
            )

            # Accumulate solution
            for M in range(NCONEC):
                JK = KCONEC_act[M, IE]
                DSOL[k] += GW[M] * DFI_act[M, IE] - HW[M] * FI_act[JK]

    return DSOL


