import numpy as np

def PROCESS_CONNECTIVITY_LINEAR(POS, QUADS):

    # POS   : (N,3) array of node coordinates
    # QUADS : (NE,10) array (9 node connectivity + physical tag)

    # Returns
    # N      : number of nodes
    # NE     : number of elements
    # KCONEC : (10, NE) array

    # Outputs
    N  = POS.shape[0]
    NE = QUADS.shape[0]

    # Slice first 9 columns, then transpose to get 9 x NE
    KCONEC = QUADS.T
    
    return N, NE, KCONEC
