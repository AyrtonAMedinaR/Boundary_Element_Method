import numpy as np

def PROCESS_CONNECTIVITY(POS, QUADS):
    # """
    # POS   : (N,3) array of node coordinates
    # QUADS : (NE,10) array (9 node connectivity + physical tag)

    # Returns
    # -------
    # N      : number of nodes
    # NE     : number of elements
    # KCONEC : (10, NE) array
    # """

    CONE = QUADS.copy()

    VECTOR_ORDER_NODES_NEW = np.array([1,5,2,6,3,7,4,8,9]) - 1
    VECTOR_ORDER_NODES_ORI = np.array([1,3,5,7,2,4,6,8,9]) - 1

    # Allocate connectivity
    Conectivity = np.zeros((CONE.shape[0], 9), dtype=int)

    # Reorder nodes
    Conectivity[:, VECTOR_ORDER_NODES_NEW] = CONE[:, VECTOR_ORDER_NODES_ORI]

    # Outputs
    N  = POS.shape[0]
    NE = Conectivity.shape[0]

    # Slice first 9 columns, then transpose to get 9 x NE
    KCONEC = Conectivity.T
    KCONEC = np.vstack([KCONEC, QUADS[:,9]])
    
    return N, NE, KCONEC
