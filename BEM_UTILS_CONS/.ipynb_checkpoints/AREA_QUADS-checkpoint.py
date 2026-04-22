import numpy as np

def AREA_QUADS(CONEC, POS, NE):

    AREA = np.zeros(NE)

    for e in range(NE):
        # node indices of element e
        n1 = CONEC[0,e]
        n2 = CONEC[1,e]
        n3 = CONEC[2,e]
        n4 = CONEC[3,e]

        # node coordinates
        P1 = POS[n1]
        P2 = POS[n2]
        P3 = POS[n3]
        P4 = POS[n4]

        # split quad into two triangles
        v1 = P2 - P1
        v2 = P3 - P1
        A1 = 0.5 * np.linalg.norm(np.cross(v1, v2))

        v3 = P4 - P1
        v4 = P3 - P1
        A2 = 0.5 * np.linalg.norm(np.cross(v3, v4))

        AREA[e] = A1 + A2

    return AREA