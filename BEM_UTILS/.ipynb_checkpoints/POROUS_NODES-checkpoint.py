import numpy as np

def POROUS_NODES(NCONEC, KCONEC_1, KCONEC_2, POS_1, POS_2, NODES_MATCH_R1R2):

    MAX_NODOS = len(NODES_MATCH_R1R2)*NCONEC
    
    STORE_NODE_R1 = np.zeros((MAX_NODOS, 3), dtype=int)
    STORE_NODE_R2 = np.zeros((MAX_NODOS, 3), dtype=int)

    REDON = 3
    L1 = 0

    # MATLAB: for J1 = NODES_MATCH_R1R2(:,3)'
    # Python: column index 2
    for J1 in NODES_MATCH_R1R2[:, 2]:

        for K1 in range(NCONEC):
            L1 += 1
            NODE = KCONEC_1[K1, J1]

            # MATLAB: for J2 = NODES_MATCH_R1R2(:,4)'
            for J2 in NODES_MATCH_R1R2[:, 3]:

                for K2 in range(NCONEC):

                    NODE2 = KCONEC_2[K2, J2]

                    if (
                        round(POS_1[NODE, 0], REDON) == round(POS_2[NODE2, 0], REDON)
                        and round(POS_1[NODE, 1], REDON) == round(POS_2[NODE2, 1], REDON)
                        and round(POS_1[NODE, 2], REDON) == round(POS_2[NODE2, 2], REDON)
                    ):
                        STORE_NODE_R1[L1-1, 0] = NODE
                        STORE_NODE_R1[L1-1, 1] = J1
                        STORE_NODE_R1[L1-1, 2] = K1

                        STORE_NODE_R2[L1-1, 0] = NODE2
                        STORE_NODE_R2[L1-1, 1] = J2
                        STORE_NODE_R2[L1-1, 2] = K2

    # MATLAB: STORE_NODE_R1(any(STORE_NODE_R1,2),:)
    mask1 = np.any(STORE_NODE_R1 != 0, axis=1)
    mask2 = np.any(STORE_NODE_R2 != 0, axis=1)

    STORE_NODE_R1U = STORE_NODE_R1[mask1]
    STORE_NODE_R2U = STORE_NODE_R2[mask2]

    return STORE_NODE_R1U, STORE_NODE_R2U