import numpy as np

# THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE NORMAL VECTORS,
# COORDINATES, JACOBIANS AND SHAPE FUNCTIONS FOR THE GAUSS POINTS

# FF → Shape functions at Gauss points
# XG1, XG2, XG3 → coordinates of Gauss points in 3-D
# ETA1, ETA2, ETA3 → unit normal vector
# XJA → surface Jacobian (area scaling)

def ETAXJA5(XJ1, XJ2, XJ3, NCONEC):
    
    # One–dimensional Gauss abscissas
    GI = np.array([
        0.238619186083197, -0.238619186083197,
        0.661209386466265, -0.661209386466265,
        0.932469514203152, -0.932469514203152
    ])

    # For a total of 36 points in a 6x6 grid
    NPG = 6
    NG = NPG * NPG

    # Shape functions
    FF = np.zeros((NCONEC, NG))

    # Shape functions derivatives
    FD1 = np.zeros(NCONEC)
    FD2 = np.zeros(NCONEC)    

    # Cartesian coordinates of Gauss points
    XG1 = np.zeros(NG)
    XG2 = np.zeros(NG)
    XG3 = np.zeros(NG)
    
    # Product of shape function derivatives
    DC  = np.zeros((NCONEC, NCONEC))
    DCT = np.zeros((NCONEC, NCONEC))

    # Surface Jacobian
    XJA  = np.zeros(NG)  

    # Normal vector at Gauss points
    ETA1 = np.zeros(NG)
    ETA2 = np.zeros(NG)
    ETA3 = np.zeros(NG)    

    for i in range(NPG):
        for j in range(NPG):

            # Order number of Gauss point
            K = NPG * i + j

            # Gauss points
            G1 = GI[i]
            G2 = GI[j]

            # Shape functions
            FF[0, K] = 0.25*(1-G1)*(1-G2)
            FF[1, K] = 0.25*(1+G1)*(1-G2)
            FF[2, K] = 0.25*(1+G1)*(1+G2)
            FF[3, K] = 0.25*(1-G1)*(1+G2)

            # Shape function derivatives dN/dxi
            FD1[0] = -0.25*(1 - G2)
            FD1[1] =  0.25*(1 - G2)
            FD1[2] =  0.25*(1 + G2)
            FD1[3] = -0.25*(1 + G2)

            # Shape function derivatives dN/deta
            FD2[0] = -0.25*(1 - G1)
            FD2[1] = -0.25*(1 + G1)
            FD2[2] =  0.25*(1 + G1)
            FD2[3] =  0.25*(1 - G1)            

            # Cartesian coordinates of Gauss point
            for ik in range(NCONEC):
                XG1[K] += FF[ik, K] * XJ1[ik]
                XG2[K] += FF[ik, K] * XJ2[ik]
                XG3[K] += FF[ik, K] * XJ3[ik]                 

            # Product of shape functions derivatives
            for IN in range(NCONEC):
                for JN in range(NCONEC):
                    DC[IN, JN]  = FD1[IN] * FD2[JN]
                    DCT[IN, JN] = FD2[IN] * FD1[JN]

            # Components of vector Jacobian - inner sum
            AUX = np.zeros(NCONEC)
            BUX = np.zeros(NCONEC)
            CUX = np.zeros(NCONEC)

            for JK in range(NCONEC):
                for IK in range(NCONEC):
                    AUX[JK] += XJ2[IK] * (DC[IK, JK] - DCT[IK, JK])
                    BUX[JK] += XJ3[IK] * (DC[IK, JK] - DCT[IK, JK])
                    CUX[JK] += XJ1[IK] * (DC[IK, JK] - DCT[IK, JK])

            # Components of vector Jacobian
            S1 = 0.0
            S2 = 0.0
            S3 = 0.0

            # Components of vector Jacobian - outer sum
            for M in range(NCONEC):
                S1 += AUX[M] * XJ3[M]
                S2 += BUX[M] * XJ1[M]
                S3 += CUX[M] * XJ2[M]

            # Surface Jacobian
            XJA[K] = np.sqrt(S1**2 + S2**2 + S3**2)

            # Normal vectors at Gauss points
            ETA1[K] = S1 / XJA[K]
            ETA2[K] = S2 / XJA[K]
            ETA3[K] = S3 / XJA[K]

    return ETA1, ETA2, ETA3, XJA, XG1, XG2, XG3, FF