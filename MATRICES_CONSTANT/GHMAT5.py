import numpy as np

from .LOCIN5 import LOCIN5
from .EXTIN5 import EXTIN5

def GHMAT5(POS, POS_M, CONE, NE, DEPTH):

    # Assemble H and G matrices for constant quadrilateral 3D BEM.

    # Parameters
    # POS_M : (NE,3) array, Collocation points (element centroids)
    # POS : (NN,3) array, Node coordinates
    # CONE : (4,NE) array, Element connectivity (zero-indexed)
    # NE : int, Number of elements

    # Returns
    # H, G : (NE,NE) complex arrays

    H = np.zeros((NE, NE), dtype=complex)
    G = np.zeros((NE, NE), dtype=complex)

    X = POS[:,0]
    Y = POS[:,1]
    Z = POS[:,2]

    XM = POS_M[:,0]
    YM = POS_M[:,1]
    ZM = POS_M[:,2]

    for J in range(NE):

        # ---- element node indices ----
        N1 = CONE[0, J]
        N2 = CONE[1, J]
        N3 = CONE[2, J]
        N4 = CONE[3, J]

        # ---- element normal vector ----
        A = (Y[N2] - Y[N1])*(Z[N3] - Z[N1]) - (Z[N2] - Z[N1])*(Y[N3] - Y[N1])
        B = (Z[N2] - Z[N1])*(X[N3] - X[N1]) - (X[N2] - X[N1])*(Z[N3] - Z[N1])
        C = (X[N2] - X[N1])*(Y[N3] - Y[N1]) - (Y[N2] - Y[N1])*(X[N3] - X[N1])

        R = np.sqrt(A*A + B*B + C*C)
        ETA = np.array([A/R, B/R, C/R])

        # ---- element corner coordinates ----
        CO = np.zeros((4,3))
        CO[0] = POS[N1]
        CO[1] = POS[N2]
        CO[2] = POS[N3]
        CO[3] = POS[N4]

        # ---- loop over collocation elements ----
        for I in range(NE):

            if I == J:
                H[I,J], G[I,J] = LOCIN5(
                    XM[I], YM[I], ZM[I], ETA,
                    X, Y, Z, N1, N2, N3, N4, DEPTH
                )
            else:
                H[I,J], G[I,J] = EXTIN5(
                    CO, XM[I], YM[I], ZM[I], ETA, DEPTH
                )

    return H, G