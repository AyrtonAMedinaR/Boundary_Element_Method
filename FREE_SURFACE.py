import numpy as np

def FREE_SURFACE(A, B, omega, g, NCONEC, COLUMNS_SI, KCONEC, NE):

    factor = (omega**2) / g

    # MATLAB: for J = COLUMNS_SI + NE
    for J in COLUMNS_SI + NE:

        # MATLAB: for K = 1:NCONEC
        for K in range(NCONEC):

            # MATLAB: IK = KCONEC(K,J)
            IK = KCONEC[K, J]

            # MATLAB: A(:,IK) = A(:,IK) - factor * B(:,(J-1)*NCONEC+K)
            A[:, IK] = A[:, IK] - factor * B[:, J*NCONEC + K]

    return A
