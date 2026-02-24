import numpy as np

def ETAXJA(XJ1, XJ2, XJ3, NCONEC):

    GI = np.array([
        0.238619186083197, -0.238619186083197,
        0.661209386466265, -0.661209386466265,
        0.932469514203152, -0.932469514203152
    ])

    NPG = 6
    NGP = NPG * NPG

    ETA1 = np.zeros(NGP)
    ETA2 = np.zeros(NGP)
    ETA3 = np.zeros(NGP)
    XJA  = np.zeros(NGP)

    XG1 = np.zeros(NGP)
    XG2 = np.zeros(NGP)
    XG3 = np.zeros(NGP)

    FF = np.zeros((9, NGP))

    K = 0
    for G1 in GI:
        for G2 in GI:

            # ---- Shape functions ----
            FF[:,K] = [
                0.25*G1*G2*(G1-1)*(G2-1),
                0.50*(1-G1**2)*G2*(G2-1),
                0.25*G1*G2*(G1+1)*(G2-1),
                0.50*G1*(G1+1)*(1-G2**2),
                0.25*G1*G2*(G1+1)*(G2+1),
                0.50*(1-G1**2)*G2*(G2+1),
                0.25*G1*G2*(G1-1)*(G2+1),
                0.50*G1*(G1-1)*(1-G2**2),
                (1-G1**2)*(1-G2**2)
            ]

            # ---- Derivatives ----
            FD1 = np.array([
                0.25*G2*(G2-1)*(2*G1-1),
                -G1*G2*(G2-1),
                0.25*G2*(G2-1)*(2*G1+1),
                0.50*(1-G2**2)*(2*G1+1),
                0.25*G2*(G2+1)*(2*G1+1),
                -G1*G2*(G2+1),
                0.25*G2*(G2+1)*(2*G1-1),
                0.50*(1-G2**2)*(2*G1-1),
                -2*G1*(1-G2**2)
            ])

            FD2 = np.array([
                0.25*G1*(G1-1)*(2*G2-1),
                0.50*(1-G1**2)*(2*G2-1),
                0.25*G1*(G1+1)*(2*G2-1),
                -G1*G2*(G1+1),
                0.25*G1*(G1+1)*(2*G2+1),
                0.50*(1-G1**2)*(2*G2+1),
                0.25*G1*(G1-1)*(2*G2+1),
                -G1*G2*(G1-1),
                -2*G2*(1-G1**2)
            ])

            # ---- Gauss coordinates ----
            XG1[K] = FF[:,K] @ XJ1
            XG2[K] = FF[:,K] @ XJ2
            XG3[K] = FF[:,K] @ XJ3

            # ---- Jacobian & normal (vectorized) ----
            a1 = XJ2 @ FD1
            a2 = XJ2 @ FD2
            b1 = XJ3 @ FD1
            b2 = XJ3 @ FD2
            c1 = XJ1 @ FD1
            c2 = XJ1 @ FD2

            AUX = a1*FD2 - a2*FD1
            BUX = b1*FD2 - b2*FD1
            CUX = c1*FD2 - c2*FD1

            S1 = AUX @ XJ3
            S2 = BUX @ XJ1
            S3 = CUX @ XJ2

            XJA[K] = np.sqrt(S1*S1 + S2*S2 + S3*S3)

            ETA1[K] = S1 / XJA[K]
            ETA2[K] = S2 / XJA[K]
            ETA3[K] = S3 / XJA[K]

            K += 1

    return ETA1, ETA2, ETA3, XJA, XG1, XG2, XG3, FF
