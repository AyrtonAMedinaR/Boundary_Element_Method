import numpy as np

def COS_DIR_LINEAR(KCONEC, J, X, Y, Z):
    # """
    # KCONEC : (10, NE) connectivity array
    # J      : element index (0-based in Python!)
    # X,Y,Z  : node coordinate arrays (size Nnodes)

    # Returns
    # -------
    # ETA : (3,) unit normal vector
    # """

    # Node numbers (convert MATLAB->Python indexing)
    N1 = KCONEC[0, J]
    N2 = KCONEC[1, J]
    N3 = KCONEC[2, J]

    A = (Y[N2] - Y[N1])*(Z[N3] - Z[N1]) - (Z[N2] - Z[N1])*(Y[N3] - Y[N1])
    B = (Z[N2] - Z[N1])*(X[N3] - X[N1]) - (X[N2] - X[N1])*(Z[N3] - Z[N1])
    C = (X[N2] - X[N1])*(Y[N3] - Y[N1]) - (Y[N2] - Y[N1])*(X[N3] - X[N1])

    R = np.sqrt(A*A + B*B + C*C)

    ETA = np.array([A/R, B/R, C/R])

    return ETA
