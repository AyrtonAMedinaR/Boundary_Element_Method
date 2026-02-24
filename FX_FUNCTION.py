import numpy as np

def FX_FUNCTION(NE_total, N_total, NCONEC, B, BC_SCAT):

    FX = np.zeros(N_total, dtype=B.dtype)

    for IEJ in range(NE_total):

        cols = IEJ * NCONEC + np.arange(NCONEC)
        FX += B[:, cols] @ BC_SCAT[:, IEJ]

    return FX
