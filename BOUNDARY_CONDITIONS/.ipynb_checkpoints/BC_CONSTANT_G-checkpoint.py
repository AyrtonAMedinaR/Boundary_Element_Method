import numpy as np

def BC_CONSTANT_G(A, B, k, NCONEC, R0_STRUCT, R1_STRUCT, N, NE, Gporous, Region_M, Region_N):

    # A and B: Coefficient matrices 
    # k: Wavenumber 
    # NCONEC: Number of nodes in a quad element 
    # R0_STRUCT: Nodes, elements and position (from 1 to 9) of structure in region 0
    # R1_STRUCT: Nodes, elements and position (from 1 to 9) of structure in region 1 
    # N: Total number of nodes
    # NE: Total number of elements
    # Gporous: Porosity parameter
    # Region_M: Outer region where the BC is applied   
    # Region_N: Inner region where the BC is applied   

    N_prev_M  = sum(N[0:Region_M])
    N_prev_N  = sum(N[0:Region_N])

    NE_prev_M = sum(NE[0:Region_M])
    NE_prev_N = sum(NE[0:Region_N])    

    # ---------- First loop ----------
    for KK in range(R0_STRUCT.shape[0]):

        NODE_R1 = R0_STRUCT[KK, 0] + N_prev_M
        ELEM_R1 = R0_STRUCT[KK, 1] + NE_prev_M
        POSI_R1 = R0_STRUCT[KK, 2]

        NODE_R2 = R1_STRUCT[KK, 0] + N_prev_N
        ELEM_R2 = R1_STRUCT[KK, 1] + NE_prev_N
        POSI_R2 = R1_STRUCT[KK, 2]

        A[:, NODE_R1] -= 1j * (Gporous * k) * (
            B[:, ELEM_R1 * NCONEC + POSI_R1] -
            B[:, ELEM_R2 * NCONEC + POSI_R2]
        )

    # ---------- Second loop ----------
    for KK in range(R1_STRUCT.shape[0]):

        NODE_R1 = R0_STRUCT[KK, 0] + N_prev_M
        ELEM_R1 = R0_STRUCT[KK, 1] + NE_prev_M
        POSI_R1 = R0_STRUCT[KK, 2]

        NODE_R2 = R1_STRUCT[KK, 0] + N_prev_N
        ELEM_R2 = R1_STRUCT[KK, 1] + NE_prev_N
        POSI_R2 = R1_STRUCT[KK, 2]

        A[:, NODE_R2] -= 1j * (Gporous * k) * (
            B[:, ELEM_R2 * NCONEC + POSI_R2] -
            B[:, ELEM_R1 * NCONEC + POSI_R1]
        )

    # ---------- First loop ----------
    for KK in range(R0_STRUCT.shape[0]):

        NODE_R1 = R0_STRUCT[KK, 0] + N_prev_M
        ELEM_R1 = R0_STRUCT[KK, 1] + NE_prev_M
        POSI_R1 = R0_STRUCT[KK, 2]

        NODE_R2 = R1_STRUCT[KK, 0] + N_prev_N
        ELEM_R2 = R1_STRUCT[KK, 1] + NE_prev_N
        POSI_R2 = R1_STRUCT[KK, 2]

        B[:, ELEM_R1 * NCONEC + POSI_R1] = 0
        B[:, ELEM_R2 * NCONEC + POSI_R2] = 0
    
    return A
