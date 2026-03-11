import numpy as np

def FX_FUNCTION(NCONEC, B, BC_S, N, NE):

    # NCONEC: Number of nodes in a quad element 
    # B: Coefficient matrix
    # BC_S: Incident velocity potential contribution in all the nodes of the incident boundary
    # N: Total number of nodes in each region
    # NE: Total number of elements in each region

    N_total  = sum(N[:]) 
    NE_total = sum(NE[:]) 

    FX = np.zeros(N_total, dtype=B.dtype)

    for IEJ in range(NE_total):

        cols = IEJ * NCONEC + np.arange(NCONEC)
        FX += B[:, cols] @ BC_S[:, IEJ]

    return FX
