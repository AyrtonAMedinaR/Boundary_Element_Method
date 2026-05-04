import numpy as np

def BC_INTERFACE(A, B, NCONEC,
                        R0_STRUCT,
                        R1_STRUCT,
                        N, NE, Region_M, Region_N):

    # A and B: Coefficient matrices 
    # NCONEC: Number of nodes in a quad element 
    # R0_STRUCT: Elements belonging to the interface in region M
    # R1_STRUCT: Element belonging to the interface in region N
    # N: Total number of nodes
    # NE: Total number of elements
    # Region_M: Region where the BC is applied
    # Region_N: Region where the BC is applied

    N_prev_M  = sum(N[0:Region_M])
    N_prev_N  = sum(N[0:Region_N])

    NE_prev_M = sum(NE[0:Region_M])
    NE_prev_N = sum(NE[0:Region_N])

    # ---- UNIQUE with 'stable' ----
    def unique_stable(arr):
        _, idx = np.unique(arr, return_index=True)
        return arr[np.sort(idx)]

    UNIQUE_VEC_1 = unique_stable(R0_STRUCT[:, 0])
    UNIQUE_VEC_2 = unique_stable(R1_STRUCT[:, 0])

    # ---- FIRST LOOP ----
    for jj in range(len(UNIQUE_VEC_1)):
        NODE_R1 = int(UNIQUE_VEC_1[jj] + N_prev_M)
        NODE_R2 = int(UNIQUE_VEC_2[jj] + N_prev_N)

        A[:, NODE_R1] += A[:, NODE_R2]

    # ---- ZERO OUT NODE_R2 IN A ----
    for jj in range(len(UNIQUE_VEC_2)):
        NODE_R2 = int(UNIQUE_VEC_2[jj] + N_prev_N)
        A[:, NODE_R2] = 0

    # ---- LOOP: MODIFY B (SUBTRACTION) ----
    for kk in range(R0_STRUCT.shape[0]):
        NODE_R1 = int(R0_STRUCT[kk, 0] + N_prev_M)
        ELEM_R1 = int(R0_STRUCT[kk, 1] + NE_prev_M)
        POSI_R1 = int(R0_STRUCT[kk, 2])

        NODE_R2 = int(R1_STRUCT[kk, 0] + N_prev_N)
        ELEM_R2 = int(R1_STRUCT[kk, 1] + NE_prev_N)
        POSI_R2 = int(R1_STRUCT[kk, 2])

        idx1 = (ELEM_R1) * NCONEC + (POSI_R1)
        idx2 = (ELEM_R2) * NCONEC + (POSI_R2)

        B[:, idx1] -= B[:, idx2]

    # ---- ZERO OUT CORRESPONDING B (R2) ----
    for kk in range(R1_STRUCT.shape[0]):
        ELEM_R2 = int(R1_STRUCT[kk, 1] + NE_prev_N)
        POSI_R2 = int(R1_STRUCT[kk, 2])

        idx2 = (ELEM_R2) * NCONEC + (POSI_R2)
        B[:, idx2] = 0

    # ---- MODIFY A USING B ----
    for kk in range(R1_STRUCT.shape[0]):
        NODE_R2 = int(R1_STRUCT[kk, 0] + N_prev_N)

        ELEM_R1 = int(R0_STRUCT[kk, 1] + NE_prev_M)
        POSI_R1 = int(R0_STRUCT[kk, 2])

        idx1 = (ELEM_R1) * NCONEC + (POSI_R1)

        A[:, NODE_R2] -= B[:, idx1]

    # ---- FINAL ZEROING OF B (R1) ----
    for kk in range(R0_STRUCT.shape[0]):
        ELEM_R1 = int(R0_STRUCT[kk, 1] + NE_prev_M)
        POSI_R1 = int(R0_STRUCT[kk, 2])

        idx1 = (ELEM_R1) * NCONEC + (POSI_R1)
        B[:, idx1] = 0

    return A, B