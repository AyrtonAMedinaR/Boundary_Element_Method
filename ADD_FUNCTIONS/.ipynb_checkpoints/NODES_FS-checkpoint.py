import numpy as np

def NODES_FS(NCONEC, ELEM_BOUNDARY, KCONEC, N, NE, Region):

    N_prev  = sum(N[0:Region])

    NE_prev = sum(NE[0:Region])
    NE_act  = sum(NE[0:Region + 1]) 

    KCONEC_M = KCONEC[:, NE_prev:NE_act] - N_prev

    total = len(ELEM_BOUNDARY) * NCONEC
    NODE_BOUNDARY = np.zeros(total, dtype=int)

    KS = 0

    for J in ELEM_BOUNDARY: # + NE_prev:
        for I in range(NCONEC):
            IK = KCONEC_M[I, J]
            NODE_BOUNDARY[KS] = IK
            KS += 1

    return NODE_BOUNDARY
