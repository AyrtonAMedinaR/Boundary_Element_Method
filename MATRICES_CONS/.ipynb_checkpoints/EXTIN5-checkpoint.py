import numpy as np
from numba import njit

from .FUNDA5 import FUNDA5

@njit(fastmath=True)
def EXTIN5(CO, XP, YP, ZP, ETA, DEPTH):

    # Compute H and G influence coefficients for a quadrilateral element
    # using 6x6 Gauss quadrature (constant element 3D BEM).

    # Parameters
    # CO : (4,3) array, Element node coordinates
    # XP, YP, ZP : float, Collocation point
    # ETA : (3,) array, Outward normal vector

    # Returns
    # H, G : complex, Influence coefficients

    OME = np.array([0.1713245, 0.3607616, 0.4679139,
                    0.4679139, 0.3607616, 0.1713245])
    GI  = np.array([-0.9324695,-0.6612094,-0.2386192,
                     0.2386192, 0.6612094, 0.9324695])

    H = 0.0 + 0j
    G = 0.0 + 0j

    # ---- outer Gauss loop ----
    for JG in range(6):
        G2 = GI[JG]
        P2 = OME[JG]

        SP = 1.0 + G2
        SM = 1.0 - G2

        # shape derivatives wrt eta
        P = np.zeros((2,4))
        P[0,0] = -0.25 * SM
        P[0,1] =  0.25 * SM
        P[0,2] =  0.25 * SP
        P[0,3] = -0.25 * SP

        # ---- inner Gauss loop ----
        for IG in range(6):
            G1 = GI[IG]
            P1 = OME[IG]

            RP = 1.0 + G1
            RM = 1.0 - G1

            # shape functions
            F = np.zeros(4)
            F[0] = 0.25 * RM * SM
            F[1] = 0.25 * RP * SM
            F[2] = 0.25 * RP * SP
            F[3] = 0.25 * RM * SP

            # shape derivatives wrt xi
            P[1,0] = -0.25 * RM
            P[1,1] = -0.25 * RP
            P[1,2] =  0.25 * RP
            P[1,3] =  0.25 * RM

            # ---- Jacobian matrix XJ (2x3) ----
            XJ = np.zeros((2,3))
            for i in range(2):
                for j in range(3):
                    for k in range(4):
                        XJ[i,j] += P[i,k] * CO[k,j]

            # ---- Jacobian determinant ----
            DET = np.sqrt(
                (XJ[0,1]*XJ[1,2] - XJ[1,1]*XJ[0,2])**2 +
                (XJ[1,0]*XJ[0,2] - XJ[0,0]*XJ[1,2])**2 +
                (XJ[0,0]*XJ[1,1] - XJ[1,0]*XJ[0,1])**2
            )

            COc = np.ascontiguousarray(CO)   # ensure contiguous once
            FG  = np.ascontiguousarray(F)
            
            XG, YG, ZG = COc.T @ FG
            
            # ---- coordinates of Gauss point ----
            #XG = np.dot(CO[:,0], F)
            #YG = np.dot(CO[:,1], F)
            #ZG = np.dot(CO[:,2], F)

            # fundamental solution
            U, Q = FUNDA5(XP, YP, ZP, ETA, XG, YG, ZG, DEPTH)

            weight = P1 * P2 * DET

            G += U * weight
            H += Q * weight

    return H, G