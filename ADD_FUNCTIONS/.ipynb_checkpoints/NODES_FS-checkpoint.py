import numpy as np

def NODES_FS(NCONEC, KCONEC, COLUMNS_R1_SI, NE0):

    total = len(COLUMNS_R1_SI) * NCONEC
    NODOS_R1_SI = np.zeros(total, dtype=int)

    KS = 0

    # MATLAB: for J = COLUMNS_R1_SI + NE0
    for J in COLUMNS_R1_SI + NE0:
        for I in range(NCONEC):
            IK = KCONEC[I, J]
            NODOS_R1_SI[KS] = IK
            KS += 1

    return NODOS_R1_SI
