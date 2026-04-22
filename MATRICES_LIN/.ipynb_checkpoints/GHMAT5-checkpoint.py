import numpy as np

from .ETAXJA5 import ETAXJA5
from .LOCIN5 import LOCIN5
from .EXTIN5 import EXTIN5

def GHMAT5(POS, NCONEC, KCONEC, NE, N, DEPTH):

    X = POS[:, 0]
    Y = POS[:, 1]
    Z = POS[:, 2]

    # Global matrices
    FI = np.zeros((N, NCONEC * NE), dtype=np.complex128)
    A  = np.zeros((N, N), dtype=np.complex128)

    # Element node coordinates
    XJ1 = np.zeros(NCONEC)
    XJ2 = np.zeros(NCONEC)
    XJ3 = np.zeros(NCONEC)

    # =====================================================
    # LOOP OVER ELEMENTS
    # =====================================================
    for IEJ in range(NE):

        # --- store element node coordinates ---
        for M in range(NCONEC):
            K = KCONEC[M, IEJ]
            XJ1[M] = X[K]
            XJ2[M] = Y[K]
            XJ3[M] = Z[K]

        # Gauss geometry for regular integration
        ETA1, ETA2, ETA3, XJA, XG1, XG2, XG3, FF = \
            ETAXJA5(XJ1, XJ2, XJ3, NCONEC)

        # =====================================================
        # LOOP OVER COLLOCATION POINTS
        # =====================================================
        for IN in range(N):

            XI = np.array([X[IN], Y[IN], Z[IN]])

            # Check if collocation node belongs to element
            NODO = -1
            for I in range(NCONEC):
                if IN == KCONEC[I, IEJ]:
                    NODO = I
                    break

            # -------------------------------------------------
            # COMPUTE ELEMENT INTEGRALS
            # -------------------------------------------------
            if NODO == -1:
                # regular integration
                HW, GW = EXTIN5(
                    ETA1, ETA2, ETA3,
                    XJA, XG1, XG2, XG3,
                    FF, XI, NCONEC, DEPTH
                )
            else:
                # singular integration (principal value)
                HW, GW = LOCIN5(
                    XJ1, XJ2, XJ3,
                    NODO, NCONEC, DEPTH
                )

            # -------------------------------------------------
            # ASSEMBLY
            # -------------------------------------------------
            for K in range(NCONEC):

                IK = KCONEC[K, IEJ]
                col = K + NCONEC * IEJ

                # ---- FI matrix ----
                FI[IN, col] += GW[K]

                # ---- A matrix ----
                if IK != IN:
                    # Off-diagonal term
                    A[IN, IK] += HW[K]

                    # Diagonal from row-sum property
                    A[IN, IN] -= HW[K]

    return A, FI