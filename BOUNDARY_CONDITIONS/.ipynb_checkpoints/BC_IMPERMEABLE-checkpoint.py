import numpy as np

def BC_IMPERMEABLE(B, NCONEC, ELEM_FS, KCONEC, NE, Region):

    # B: Coefficient matrices 
    # NCONEC: Number of nodes in a quad element 
    # ELEM_FS: Elements belonging to the free surface
    # KCONEC: Element connectivity array
    # NE: Total number of elements in each region
    # Region: Region where the BC is applied

    NE_prev = sum(NE[0:Region])

    for J in ELEM_FS + NE_prev:
        for K in range(NCONEC):

            B[:, J*NCONEC + K] = 0 
  
    return B
