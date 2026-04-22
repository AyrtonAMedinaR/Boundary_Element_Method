import numpy as np

def PROCESS_CONNECTIVITY_CONSTANT(MESH_POS, CONE):
    
    # Compute element centroids and connectivity in column format.

    # Parameters
    # CONE : (NE,4) array, Element connectivity (zero-indexed)
    # MESH_POS : (N,3) array, Node coordinates

    # Returns
    # MESH_POS_M : (NE,3) array, Element centroids
    # KCONEC : (4,NE) array, Transposed connectivity
    # N : int, Number of nodes
    # NE : int, Number of elements

    N  = MESH_POS.shape[0]
    NE = CONE.shape[0]

    X = MESH_POS[:,0]
    Y = MESH_POS[:,1]
    Z = MESH_POS[:,2]

    XM = np.zeros(NE)
    YM = np.zeros(NE)
    ZM = np.zeros(NE)

    for i in range(NE):

        N1 = CONE[i,0]
        N2 = CONE[i,1]
        N3 = CONE[i,2]
        N4 = CONE[i,3]

        # Quadrilateral centroid
        XM[i] = (X[N1] + X[N2] + X[N3] + X[N4]) / 4.0
        YM[i] = (Y[N1] + Y[N2] + Y[N3] + Y[N4]) / 4.0
        ZM[i] = (Z[N1] + Z[N2] + Z[N3] + Z[N4]) / 4.0

    MESH_POS_M = np.column_stack((XM, YM, ZM))

    # Transpose connectivity to match BEM assembly routine
    KCONEC = CONE.T

    return N, NE, KCONEC, MESH_POS_M