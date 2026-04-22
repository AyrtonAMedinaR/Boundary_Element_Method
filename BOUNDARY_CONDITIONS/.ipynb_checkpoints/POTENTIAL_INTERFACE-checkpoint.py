import numpy as np

def POTENTIAL_INTERFACE(PHI_ONLY, X_Vec, R0_STRUCT, R1_STRUCT, N, Region_M, Region_N):

    N_prev_M  = sum(N[0:Region_M])
    N_prev_N  = sum(N[0:Region_N])
    
    for jj in range(len(R0_STRUCT)):
        NODE_R0 = int(R0_STRUCT[jj,0] + N_prev_M)
        NODE_R1 = int(R1_STRUCT[jj,0] + N_prev_N)
    
        PHI_ONLY[NODE_R1] = X_Vec[NODE_R0]

    return PHI_ONLY
