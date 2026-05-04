import numpy as np

def POTENTIAL_INTERFACE(PHI_ONLY, Coef, R0_STRUCT, R1_STRUCT, N, Region_M, Region_N):

    # PHI_ONLY: Vector containing only PHI values
    # Coef: Vector containing PHI and DPHI values
    # R0_STRUCT: Elements belonging to the interface in region M
    # R1_STRUCT: Element belonging to the interface in region N
    # N: Total number of nodes
    # Region_M: Region where the BC is applied
    # Region_N: Region where the BC is applied    

    N_prev_M  = sum(N[0:Region_M])
    N_prev_N  = sum(N[0:Region_N])
    
    for jj in range(len(R0_STRUCT)):
        NODE_R0 = int(R0_STRUCT[jj,0]) + N_prev_M
        NODE_R1 = int(R1_STRUCT[jj,0]) + N_prev_N
    
        PHI_ONLY[NODE_R1] = Coef[NODE_R0]

    return PHI_ONLY
