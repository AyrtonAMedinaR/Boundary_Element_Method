import numpy as np

def COORDINATES(POS0, QUADS0):
    # """
    # POS0   : (N0,3) array of node coordinates
    # QUADS0 : (NE0,10) array (9 node connectivity + physical tag)

    # Returns
    # -------
    # N0      : number of nodes
    # NE0     : number of elements
    # KCONEC0 : (10, NE0) array
    # """

    CONE0 = QUADS0.copy()

    VECTOR_ORDER_NODES_NEW = np.array([1,5,2,6,3,7,4,8,9]) - 1
    VECTOR_ORDER_NODES_ORI = np.array([1,3,5,7,2,4,6,8,9]) - 1

    # Allocate connectivity
    Conectivity0 = np.zeros((CONE0.shape[0], 9), dtype=int)

    # Reorder nodes
    Conectivity0[:, VECTOR_ORDER_NODES_NEW] = CONE0[:, VECTOR_ORDER_NODES_ORI]

    # Outputs
    N0  = POS0.shape[0]
    NE0 = Conectivity0.shape[0]

    KCONEC0 = Conectivity0.T
    KCONEC0 = np.vstack([KCONEC0, QUADS0[:,9]])

    return N0, NE0, KCONEC0
