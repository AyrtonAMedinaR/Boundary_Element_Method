import numpy as np

def FAR_FIELD(A, B, k, NCONEC, COLUMNS_FF, KCONEC, NE):

    for J in COLUMNS_FF + NE:

        IKs = KCONEC[:, J]                 # shape (NCONEC,)
        cols_B = J*NCONEC + np.arange(NCONEC)

        A[:, IKs] -= 1j * k * B[:, cols_B]

    return A

