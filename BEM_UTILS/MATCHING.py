import numpy as np

def MATCHING(REF_1, POS_1, KCONEC_1,
             REF_2, POS_2, KCONEC_2):

    REDON = 3

    # ---- BUILD BNL_1 ----
    n1 = POS_1.shape[0]
    BNL_1 = np.zeros((n1,4))
    BNL_1[:,0] = np.arange(0, n1)
    BNL_1[:,1] = np.round(POS_1[:,0], REDON)
    BNL_1[:,2] = np.round(POS_1[:,1], REDON)
    BNL_1[:,3] = np.round(POS_1[:,2], REDON)

    # ---- BUILD BNL_2 ----
    n2 = POS_2.shape[0]
    BNL_2 = np.zeros((n2,4))
    BNL_2[:,0] = np.arange(0, n2)
    BNL_2[:,1] = np.round(POS_2[:,0], REDON)
    BNL_2[:,2] = np.round(POS_2[:,1], REDON)
    BNL_2[:,3] = np.round(POS_2[:,2], REDON)

    # ---- NODES INVOLVED IN INTERFACE ----
    POS_1_R1R2 = KCONEC_1[8, REF_1]
    POS_2_R2R1 = KCONEC_2[8, REF_2]

    NODES_MATCH = np.zeros((len(POS_1_R1R2),4), dtype=int)

    JJJ = 0

    for J in POS_1_R1R2:
        for K in POS_2_R2R1:
            if np.all(BNL_1[J,1:4] == BNL_2[K,1:4]):
                NODES_MATCH[JJJ,0] = BNL_1[J,0]
                NODES_MATCH[JJJ,1] = BNL_2[K,0]
                JJJ += 1

    # ---- IDENTIFY CONNECTED ELEMENTS ----

    for J in range(int(np.min(REF_1)), int(np.max(REF_1))+1):
        for K in range(len(POS_1_R1R2)):
            if KCONEC_1[8,J] == NODES_MATCH[K,0]:
                NODES_MATCH[K,2] = J

    for J in range(int(np.min(REF_2)), int(np.max(REF_2))+1):
        for K in range(len(POS_1_R1R2)):
            if KCONEC_2[8,J] == NODES_MATCH[K,1]:
                NODES_MATCH[K,3] = J

    return NODES_MATCH
