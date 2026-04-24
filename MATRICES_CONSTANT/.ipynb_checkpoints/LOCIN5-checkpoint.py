import numpy as np
from numba import njit

from .EXTIN5 import EXTIN5

@njit(fastmath=True)
def LOCIN5(XM, YM, ZM, ETA, X, Y, Z, N1, N2, N3, N4, DEPTH):

    # Computes diagonal coefficient of G matrix (singular integral)
    # using collapsed quadrilateral -> 4 internal triangles.

    # Parameters
    # XM, YM, ZM : float, Collocation point (element centroid)
    # ETA : (3,) array, Element normal
    # X, Y, Z : arrays, Global node coordinates (zero-indexed)
    # N1, N2, N3, N4 : int, Element node indices (zero-indexed)

    # Returns
    # H, G : complex

    CO = np.zeros((4,3))

    # Nodes 3 and 4 collapse to collocation point
    CO[2] = [XM, YM, ZM]
    CO[3] = [XM, YM, ZM]

    G = 0.0 + 0j

    for IT in range(4):

        if IT == 0:
            CO[0] = [X[N1], Y[N1], Z[N1]]
            CO[1] = [X[N2], Y[N2], Z[N2]]

        elif IT == 1:
            CO[0] = [X[N2], Y[N2], Z[N2]]
            CO[1] = [X[N3], Y[N3], Z[N3]]

        elif IT == 2:
            CO[0] = [X[N3], Y[N3], Z[N3]]

            if N3 != N4:
                CO[1] = [X[N4], Y[N4], Z[N4]]
            else:
                CO[1] = [X[N1], Y[N1], Z[N1]]

        elif IT == 3:
            if N3 != N4:
                CO[0] = [X[N4], Y[N4], Z[N4]]
                CO[1] = [X[N1], Y[N1], Z[N1]]
            else:
                continue  # matches MATLAB behavior (no integration)

        # integrate collapsed quadrilateral
        _, GP = EXTIN5(CO, XM, YM, ZM, ETA, DEPTH)
        G += GP

    H = 0.5 + 0j

    return H, G