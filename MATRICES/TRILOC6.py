import numpy as np
from numba import njit

from .FUNDA6 import FUNDA6

@njit(fastmath=True)
def TRILOC6(HW, GW, XI, XJ1, XJ2, XJ3, AREA, V1, V2, NCONEC, DEPTH):

    # One dimensional Gauss abscissas
    T = np.array([
        0.238619186083197, -0.238619186083197,
        0.661209386466265, -0.661209386466265,
        0.932469514203152, -0.932469514203152
    ])

    # Gauss weights
    OME = np.array([
        0.467913934572691, 0.467913934572691,
        0.360761573048139, 0.360761573048139,
        0.171324492379170, 0.171324492379170
    ])

    # Shape functions
    F   = np.zeros(9)

    # Shape functions derivatives
    FD1 = np.zeros(9)
    FD2 = np.zeros(9)

    # Terms belonging to the components of the Jacobian determinant
    AUX = np.zeros(NCONEC)
    BUX = np.zeros(NCONEC)
    CUX = np.zeros(NCONEC)

    # Product of shape function derivatives
    DC  = np.zeros((NCONEC, NCONEC))
    DCT = np.zeros((NCONEC, NCONEC))

    # Normal vector at Gauss points
    ETA = np.zeros(3)

    # Cartesian coordinates of Gauss points
    X   = np.zeros(3)

    # Product of shape function, surface Jacobian and weight
    P12 = np.zeros(NCONEC)

    for i in range(6):
        for j in range(6):

            SS1 = 0.5 * (T[i] + 1.0)
            SS2 = 0.5 * (T[j] + 1.0)

            G1 = (1.0 - SS1)*V1[0] + SS1*(1.0 - SS2)*V1[1] + SS1*SS2*V1[2]
            G2 = (1.0 - SS1)*V2[0] + SS1*(1.0 - SS2)*V2[1] + SS1*SS2*V2[2]

            # Shape functions
            F[0] = 0.25*G1*G2*(G1-1.0)*(G2-1.0)
            F[1] = 0.50*(1-G1*G1)*G2*(G2-1.0)
            F[2] = 0.25*G1*G2*(G1+1.0)*(G2-1.0)
            F[3] = 0.50*G1*(G1+1.0)*(1-G2*G2)
            F[4] = 0.25*G1*G2*(G1+1.0)*(G2+1.0)
            F[5] = 0.50*(1-G1*G1)*G2*(G2+1.0)
            F[6] = 0.25*G1*G2*(G1-1.0)*(G2+1.0)
            F[7] = 0.50*G1*(G1-1.0)*(1-G2*G2)
            F[8] = (1-G1*G1)*(1-G2*G2)

            # Shape functions derivatives
            FD1[0] = 0.25*G2*(G2-1.0)*(2*G1-1.0)
            FD1[1] = -G1*G2*(G2-1.0)
            FD1[2] = 0.25*G2*(G2-1.0)*(2*G1+1.0)
            FD1[3] = 0.50*(1-G2*G2)*(2*G1+1.0)
            FD1[4] = 0.25*G2*(G2+1.0)*(2*G1+1.0)
            FD1[5] = -G1*G2*(G2+1.0)
            FD1[6] = 0.25*G2*(G2+1.0)*(2*G1-1.0)
            FD1[7] = 0.50*(1-G2*G2)*(2*G1-1.0)
            FD1[8] = -2.0*G1*(1-G2*G2)

            FD2[0] = 0.25*G1*(G1-1.0)*(2*G2-1.0)
            FD2[1] = 0.50*(1-G1*G1)*(2*G2-1.0)
            FD2[2] = 0.25*G1*(G1+1.0)*(2*G2-1.0)
            FD2[3] = -G1*G2*(G1+1.0)
            FD2[4] = 0.25*G1*(G1+1.0)*(2*G2+1.0)
            FD2[5] = 0.50*(1-G1*G1)*(2*G2+1.0)
            FD2[6] = 0.25*G1*(G1-1.0)*(2*G2+1.0)
            FD2[7] = -G1*G2*(G1-1.0)
            FD2[8] = -2.0*G2*(1-G1*G1)

            # Reset of terms of the Jacobian determinant
            for k in range(NCONEC):
                AUX[k] = 0.0
                BUX[k] = 0.0
                CUX[k] = 0.0

            # Product of shape function derivatives
            for IN in range(NCONEC):
                for JN in range(NCONEC):
                    DC[IN,JN]  = FD1[IN]*FD2[JN]
                    DCT[IN,JN] = FD2[IN]*FD1[JN]

            for JK in range(NCONEC):
                for IK in range(NCONEC):
                    diff = DC[IK,JK] - DCT[IK,JK]
                    AUX[JK] += XJ2[IK]*diff
                    BUX[JK] += XJ3[IK]*diff
                    CUX[JK] += XJ1[IK]*diff

            S1 = 0.0
            S2 = 0.0
            S3 = 0.0

            # Components of the Jacobian determinant
            for M in range(NCONEC):
                S1 += AUX[M]*XJ3[M]
                S2 += BUX[M]*XJ1[M]
                S3 += CUX[M]*XJ2[M]

            # Surface Jacobian
            XJA = np.sqrt(S1*S1 + S2*S2 + S3*S3)

            # Normal vector at Gauss points
            ETA[0] = S1 / XJA
            ETA[1] = S2 / XJA
            ETA[2] = S3 / XJA

            X[0] = 0.0
            X[1] = 0.0
            X[2] = 0.0

            # Cartesian coordinates of Gauss points            
            for IK in range(NCONEC):
                X[0] += F[IK]*XJ1[IK]
                X[1] += F[IK]*XJ2[IK]
                X[2] += F[IK]*XJ3[IK]

            # Fundamental solution
            U, Q = FUNDA6(XI, ETA, X, DEPTH)

            # Weight
            weight = (T[i]+1.0)*AREA*0.25*OME[i]*OME[j]

            for IK in range(NCONEC):
                P12[IK] = F[IK]*XJA*weight

            # Summation for integrals
            for KK in range(NCONEC):
                GW[KK] += U*P12[KK]
                HW[KK] += Q*P12[KK]

    return HW, GW