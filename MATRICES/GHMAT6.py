import numpy as np

from .ETAXJA import ETAXJA
from .LOCIN6 import LOCIN6
from .EXTIN6 import EXTIN6

def GHMAT6(POS, NCONEC, KCONEC, NE, N, hd):

    X = POS[:, 0]
    Y = POS[:, 1]
    Z = POS[:, 2]

    FI = np.zeros((N, NCONEC * NE), dtype=np.complex128)
    A  = np.zeros((N, N), dtype=np.complex128)

    HS11 = 0.0 + 0.0j

    XJ1 = np.zeros(NCONEC)
    XJ2 = np.zeros(NCONEC)
    XJ3 = np.zeros(NCONEC)

    # --------------------------------------------------
    for IEJ in range(NE):

        # Store element nodes
        for M in range(NCONEC):
            K = KCONEC[M, IEJ]
            XJ1[M] = X[K]
            XJ2[M] = Y[K]
            XJ3[M] = Z[K]

        ETA1, ETA2, ETA3, XJA, XG1, XG2, XG3, FF = \
            ETAXJA(XJ1, XJ2, XJ3, NCONEC)

        # --------------------------------------------------
        for IN in range(N):

            XI = np.array([X[IN], Y[IN], Z[IN]])

            # ---- Python-style NODO ----
            NODO = -1
            for I in range(NCONEC):
                if IN == KCONEC[I, IEJ]:
                    NODO = I
                    break

            # ---- Integrals ----
            if NODO == -1:
                HW, GW, HS = EXTIN6(
                    ETA1, ETA2, ETA3,
                    XJA, XG1, XG2, XG3,
                    FF, XI, NCONEC, hd
                )

            else:
                HW, GW, HS, HDIF = LOCIN6(
                    XJ1, XJ2, XJ3,
                    NODO, NCONEC, hd
                )

            # ---- Diagonal terms ----
            D1 = 0.0 + 0.0j
            D2 = 0.0 + 0.0j

            if NODO == -1:
                for I in range(NCONEC):
                    D1 -= HS[I]
                if IN == 0:
                    HS11 += D1
                    
            else:
                for I in range(NCONEC):
                    if I == NODO:
                        D2 += HDIF
                    else:
                        D2 -= HS[I]
                        if IN == 0:
                            HS11 -= HS[I]

            # ---- Assembly ----
            for K in range(NCONEC):

                IK = KCONEC[K, IEJ]
                col = K + NCONEC * IEJ
                
                FI[IN, col] += GW[K]         

                if IN == IK:
                    A[IN, IN] += D2
                else:
                    A[IN, IK] += HW[K]

            # ---- Extra diagonal ----
            if NODO == -1:
                A[IN, IN] += D1

        # end IN
    # end IEJ

    return A, FI
