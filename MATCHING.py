import numpy as np

def MATCHING(COLUMNS_R1_R2, POS1, KCONEC1,
             COLUMNS_R2_R1, POS2, KCONEC2):

    REDON = 3

    # ---- BUILD BNL1 ----
    n1 = POS1.shape[0]
    BNL1 = np.zeros((n1,4))
    BNL1[:,0] = np.arange(0, n1)
    BNL1[:,1] = np.round(POS1[:,0], REDON)
    BNL1[:,2] = np.round(POS1[:,1], REDON)
    BNL1[:,3] = np.round(POS1[:,2], REDON)

    # ---- BUILD BNL2 ----
    n2 = POS2.shape[0]
    BNL2 = np.zeros((n2,4))
    BNL2[:,0] = np.arange(0, n2)
    BNL2[:,1] = np.round(POS2[:,0], REDON)
    BNL2[:,2] = np.round(POS2[:,1], REDON)
    BNL2[:,3] = np.round(POS2[:,2], REDON)

    # ---- NODES INVOLVED IN INTERFACE ----
    POS1_R1R2 = KCONEC1[8, COLUMNS_R1_R2]
    POS2_R2R1 = KCONEC2[8, COLUMNS_R2_R1]

    NODES_MATCH = np.zeros((len(POS1_R1R2),4), dtype=int)

    JJJ = 0

    for J in POS1_R1R2:
        for K in POS2_R2R1:
            if np.all(BNL1[J,1:4] == BNL2[K,1:4]):
                NODES_MATCH[JJJ,0] = BNL1[J,0]
                NODES_MATCH[JJJ,1] = BNL2[K,0]
                JJJ += 1

    # ---- IDENTIFY CONNECTED ELEMENTS ----

    for J in range(int(np.min(COLUMNS_R1_R2)), int(np.max(COLUMNS_R1_R2))+1):
        for K in range(len(POS1_R1R2)):
            if KCONEC1[8,J] == NODES_MATCH[K,0]:
                NODES_MATCH[K,2] = J

    for J in range(int(np.min(COLUMNS_R2_R1)), int(np.max(COLUMNS_R2_R1))+1):
        for K in range(len(POS1_R1R2)):
            if KCONEC2[8,J] == NODES_MATCH[K,1]:
                NODES_MATCH[K,3] = J

    return NODES_MATCH
