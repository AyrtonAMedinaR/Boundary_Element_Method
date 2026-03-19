import numpy as np

# THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE NORMAL VECTORS,
# COORDINATES, JACOBIANS AND SHAPE FUNCTIONS FOR THE GAUSS POINTS

def ETAXJA(XJ1, XJ2, XJ3, NCONEC):
    
    # GAUSS POINTS FOR QUADRILATERAL ELEMENTS
    GI = np.array([
        0.238619186083197, -0.238619186083197,
        0.661209386466265, -0.661209386466265,
        0.932469514203152, -0.932469514203152
    ])

    NPG = 6
    NG = NPG * NPG

    ETA1 = np.zeros(NG)
    ETA2 = np.zeros(NG)
    ETA3 = np.zeros(NG)
    XJA  = np.zeros(NG)

    XG1 = np.zeros(NG)
    XG2 = np.zeros(NG)
    XG3 = np.zeros(NG)

    FF = np.zeros((9, NG))

    DC  = np.zeros((NCONEC, NCONEC))
    DCT = np.zeros((NCONEC, NCONEC))

    FD1 = np.zeros(9)
    FD2 = np.zeros(9)

    for i in range(NPG):
        for j in range(NPG):

            # Order number of Gauss point
            K = NPG * i + j

            G1 = GI[i]
            G2 = GI[j]

            # Shape functions
            FF[0, K] = 0.25*G1*G2*(G1-1.0)*(G2-1.0)
            FF[1, K] = 0.50*(1.0-G1**2)*G2*(G2-1.0)
            FF[2, K] = 0.25*G1*G2*(G1+1.0)*(G2-1.0)
            FF[3, K] = 0.50*G1*(G1+1.0)*(1.0-G2**2)
            FF[4, K] = 0.25*G1*G2*(G1+1.0)*(G2+1.0)
            FF[5, K] = 0.50*(1.0-G1**2)*G2*(G2+1.0)
            FF[6, K] = 0.25*G1*G2*(G1-1.0)*(G2+1.0)
            FF[7, K] = 0.50*G1*(G1-1.0)*(1.0-G2**2)
            FF[8, K] = (1.0-G1**2)*(1.0-G2**2)

            # Shape function derivatives
            FD1[0] = 0.25*G2*(G2-1.0)*(2*G1-1.0)
            FD1[1] = -G1*G2*(G2-1.0)
            FD1[2] = 0.25*G2*(G2-1.0)*(2*G1+1.0)
            FD1[3] = 0.50*(1.0-G2**2)*(2*G1+1.0)
            FD1[4] = 0.25*G2*(G2+1.0)*(2*G1+1.0)
            FD1[5] = -G1*G2*(G2+1.0)
            FD1[6] = 0.25*G2*(G2+1.0)*(2*G1-1.0)
            FD1[7] = 0.50*(1.0-G2**2)*(2*G1-1.0)
            FD1[8] = -2.0*G1*(1.0-G2**2)

            FD2[0] = 0.25*G1*(G1-1.0)*(2*G2-1.0)
            FD2[1] = 0.50*(1.0-G1**2)*(2*G2-1.0)
            FD2[2] = 0.25*G1*(G1+1.0)*(2*G2-1.0)
            FD2[3] = -G1*G2*(G1+1.0)
            FD2[4] = 0.25*G1*(G1+1.0)*(2*G2+1.0)
            FD2[5] = 0.50*(1.0-G1**2)*(2*G2+1.0)
            FD2[6] = 0.25*G1*(G1-1.0)*(2*G2+1.0)
            FD2[7] = -G1*G2*(G1-1.0)
            FD2[8] = -2.0*G2*(1.0-G1**2)

            # Cartesian coordinates of Gauss point
            for ik in range(NCONEC):
                XG1[K] += FF[ik, K] * XJ1[ik]
                XG2[K] += FF[ik, K] * XJ2[ik]
                XG3[K] += FF[ik, K] * XJ3[ik]

            # Jacobian and normal vector
            for IN in range(NCONEC):
                for JN in range(NCONEC):
                    DC[IN, JN]  = FD1[IN] * FD2[JN]
                    DCT[IN, JN] = FD2[IN] * FD1[JN]

            AUX = np.zeros(NCONEC)
            BUX = np.zeros(NCONEC)
            CUX = np.zeros(NCONEC)

            for JK in range(NCONEC):
                for IK in range(NCONEC):
                    AUX[JK] += XJ2[IK] * (DC[IK, JK] - DCT[IK, JK])
                    BUX[JK] += XJ3[IK] * (DC[IK, JK] - DCT[IK, JK])
                    CUX[JK] += XJ1[IK] * (DC[IK, JK] - DCT[IK, JK])

            S1 = 0.0
            S2 = 0.0
            S3 = 0.0

            for M in range(NCONEC):
                S1 += AUX[M] * XJ3[M]
                S2 += BUX[M] * XJ1[M]
                S3 += CUX[M] * XJ2[M]

            XJA[K] = np.sqrt(S1**2 + S2**2 + S3**2)

            ETA1[K] = S1 / XJA[K]
            ETA2[K] = S2 / XJA[K]
            ETA3[K] = S3 / XJA[K]

    return ETA1, ETA2, ETA3, XJA, XG1, XG2, XG3, FF