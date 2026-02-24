import numpy as np
from numba import njit
from FUNDA6 import FUNDA6

# IMPORTANT:
# This function assumes FUNDA6_nb is already defined with @njit

@njit(fastmath=False)
def EXTIN6(ETA1, ETA2, ETA3,
              XJA,
              XG1, XG2, XG3,
              FF,
              XI,
              NCONEC,
              hd):

    # ---- Gauss weights (hardcoded for Numba)
    OME = np.array((
        0.467913934572691, 0.467913934572691,
        0.360761573048139, 0.360761573048139,
        0.171324492379170, 0.171324492379170
    ))

    HW = np.zeros(NCONEC, dtype=np.complex128)
    GW = np.zeros(NCONEC, dtype=np.complex128)
    HS = np.zeros(NCONEC, dtype=np.complex128)  # float64

    NPG = 6

    for i in range(NPG):
        for j in range(NPG):

            K = NPG*i + j

            # --- normal components directly (avoid array creation)
            eta1 = ETA1[K]
            eta2 = ETA2[K]
            eta3 = ETA3[K]

            # ---------- weight
            w = XJA[K] * OME[i] * OME[j]

            # ---------- fundamental solution
            U, Q, QS = FUNDA6(
                XI,
                np.array((eta1, eta2, eta3)),
                XG1, XG2, XG3,
                K,
                hd
            )

            for n in range(NCONEC):
                P12 = FF[n, K] * w

                GW[n] += U  * P12
                HW[n] += Q  * P12
                HS[n] += QS * P12

    return HW, GW, HS










