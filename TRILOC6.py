import numpy as np
from numba import njit

# FUNDA6_TRI MUST also be njit compiled
from FUNDA6_TRI import FUNDA6_TRI

@njit(cache=True)
def TRILOC6(HW, GW, HS, XI, XJ1, XJ2, XJ3, AREA, V1, V2, NCONEC, hd):

    T = np.array([
        0.238619186083197, -0.238619186083197,
        0.661209386466265, -0.661209386466265,
        0.932469514203152, -0.932469514203152
    ])

    OME = np.array([
        0.467913934572691, 0.467913934572691,
        0.360761573048139, 0.360761573048139,
        0.171324492379170, 0.171324492379170
    ])

    for i in range(6):
        for j in range(6):

            SS1 = 0.5*(T[i] + 1.0)
            SS2 = 0.5*(T[j] + 1.0)

            G1 = (1.0-SS1)*V1[0] + SS1*(1.0-SS2)*V1[1] + SS1*SS2*V1[2]
            G2 = (1.0-SS1)*V2[0] + SS1*(1.0-SS2)*V2[1] + SS1*SS2*V2[2]

            # ---------------- Shape functions
            F = np.empty(9)
            FD1 = np.empty(9)
            FD2 = np.empty(9)

            F[0] = 0.25*G1*G2*(G1-1)*(G2-1)
            F[1] = 0.5*(1-G1*G1)*G2*(G2-1)
            F[2] = 0.25*G1*G2*(G1+1)*(G2-1)
            F[3] = 0.5*G1*(G1+1)*(1-G2*G2)
            F[4] = 0.25*G1*G2*(G1+1)*(G2+1)
            F[5] = 0.5*(1-G1*G1)*G2*(G2+1)
            F[6] = 0.25*G1*G2*(G1-1)*(G2+1)
            F[7] = 0.5*G1*(G1-1)*(1-G2*G2)
            F[8] = (1-G1*G1)*(1-G2*G2)

            FD1[0] = 0.25*G2*(G2-1)*(2*G1-1)
            FD1[1] = -G1*G2*(G2-1)
            FD1[2] = 0.25*G2*(G2-1)*(2*G1+1)
            FD1[3] = 0.5*(1-G2*G2)*(2*G1+1)
            FD1[4] = 0.25*G2*(G2+1)*(2*G1+1)
            FD1[5] = -G1*G2*(G2+1)
            FD1[6] = 0.25*G2*(G2+1)*(2*G1-1)
            FD1[7] = 0.5*(1-G2*G2)*(2*G1-1)
            FD1[8] = -2*G1*(1-G2*G2)

            FD2[0] = 0.25*G1*(G1-1)*(2*G2-1)
            FD2[1] = 0.5*(1-G1*G1)*(2*G2-1)
            FD2[2] = 0.25*G1*(G1+1)*(2*G2-1)
            FD2[3] = -G1*G2*(G1+1)
            FD2[4] = 0.25*G1*(G1+1)*(2*G2+1)
            FD2[5] = 0.5*(1-G1*G1)*(2*G2+1)
            FD2[6] = 0.25*G1*(G1-1)*(2*G2+1)
            FD2[7] = -G1*G2*(G1-1)
            FD2[8] = -2*G2*(1-G1*G1)

            # ---------- Dot products
            a1 = 0.0
            a2 = 0.0
            b1 = 0.0
            b2 = 0.0
            c1 = 0.0
            c2 = 0.0

            for k in range(9):
                a1 += XJ2[k]*FD1[k]
                a2 += XJ2[k]*FD2[k]
                b1 += XJ3[k]*FD1[k]
                b2 += XJ3[k]*FD2[k]
                c1 += XJ1[k]*FD1[k]
                c2 += XJ1[k]*FD2[k]

            AUX = a1*FD2 - a2*FD1
            BUX = b1*FD2 - b2*FD1
            CUX = c1*FD2 - c2*FD1

            S1 = 0.0
            S2 = 0.0
            S3 = 0.0

            for k in range(9):
                S1 += AUX[k]*XJ3[k]
                S2 += BUX[k]*XJ1[k]
                S3 += CUX[k]*XJ2[k]

            XJA = np.sqrt(S1*S1 + S2*S2 + S3*S3)

            ETA0 = S1/XJA
            ETA1 = S2/XJA
            ETA2 = S3/XJA

            X0 = 0.0
            X1 = 0.0
            X2 = 0.0

            for k in range(9):
                X0 += F[k]*XJ1[k]
                X1 += F[k]*XJ2[k]
                X2 += F[k]*XJ3[k]

            U, Q, QS = FUNDA6_TRI(
                XI[0], XI[1], XI[2],
                ETA0, ETA1, ETA2,
                X0, X1, X2, hd
            )

            weight = XJA*(T[i]+1.0)*AREA*0.25*OME[i]*OME[j]

            for k in range(9):
                val = F[k]*weight
                GW[k] += U*val
                HW[k] += Q*val
                HS[k] += QS*val

    return HW, GW, HS
