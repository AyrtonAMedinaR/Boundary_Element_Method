import numpy as np
from numba import njit

from .FUNDA5 import FUNDA5

@njit(fastmath=True)
def EXTIN5(ETA1, ETA2, ETA3,
              XJA,
              XG1, XG2, XG3,
              FF,
              XI,
              NCONEC, DEPTH):

    # Gauss weights
    OME = np.array((
        0.467913934572691, 0.467913934572691,
        0.360761573048139, 0.360761573048139,
        0.171324492379170, 0.171324492379170
    ))

    # Element integrals
    HW = np.zeros(NCONEC, dtype=np.complex128)
    GW = np.zeros(NCONEC, dtype=np.complex128)

    # For a total of 36 points in a 6x6 grid
    NPG = 6

    for i in range(NPG):
        for j in range(NPG):

            # Order number of Gauss point
            K = NPG*i + j

            # Normal components
            eta1 = ETA1[K]
            eta2 = ETA2[K]
            eta3 = ETA3[K]

            # Weight
            w = XJA[K] * OME[i] * OME[j]

            # Fundamental solution
            U, Q = FUNDA5(
                XI,
                np.array((eta1, eta2, eta3)),
                np.array((XG1[K], XG2[K], XG3[K])), 
                DEPTH
            )

            for n in range(NCONEC):
                P12 = FF[n, K] * w

                # Summation for integrals
                GW[n] += U  * P12
                HW[n] += Q  * P12

    return HW, GW










