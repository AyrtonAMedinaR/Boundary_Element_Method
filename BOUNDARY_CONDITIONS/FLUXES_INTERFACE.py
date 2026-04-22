import numpy as np

def FLUXES_INTERFACE(PHI, DPHI,
                        R0_STRUCT,
                        R1_STRUCT,
                        N, NE, Region_M, Region_N):

    N_prev_M  = sum(N[0:Region_M])
    N_prev_N  = sum(N[0:Region_N])

    NE_prev_M = sum(NE[0:Region_M])
    NE_prev_N = sum(NE[0:Region_N])
    
    for kk in range(R0_STRUCT.shape[0]):
        NODE_R0 = int(R0_STRUCT[kk, 0]) + N_prev_M
        ELEM_R0 = int(R0_STRUCT[kk, 1]) + NE_prev_M
        POSI_R0 = int(R0_STRUCT[kk, 2])

        NODE_R1 = int(R1_STRUCT[kk, 0]) + N_prev_N
        ELEM_R1 = int(R1_STRUCT[kk, 1]) + NE_prev_N
        POSI_R1 = int(R1_STRUCT[kk, 2])
        
        DPHI[POSI_R0, ELEM_R0] =   PHI[NODE_R1]
        DPHI[POSI_R1, ELEM_R1] = - PHI[NODE_R1]

    return DPHI
