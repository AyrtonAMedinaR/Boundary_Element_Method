import numpy as np

def ELEMENT_CONNECTION(NCONEC, KCONEC, POS, NE):
    # """
    # NCONEC  : number of nodes per element (9)
    # KCONEC : (10, NE) array
    # POS    : (N1,3) array
    # NE     : number of elements

    # Returns
    # -------
    # REF_ELEM : dictionary of element groups
    # NORMAL  : (NE,3) array of normals
    # """

    from .COS_DIR import COS_DIR

    NORMAL = np.zeros((NE, 3))

    # Coordinates
    X = POS[:,0]
    Y = POS[:,1]
    Z = POS[:,2]

    # ---- NORMAL VECTORS ----
    for J in range(NE):
        ETA = COS_DIR(KCONEC, J, X, Y, Z)

        NORMAL[J,0] = ETA[0]
        NORMAL[J,1] = ETA[1]
        NORMAL[J,2] = ETA[2]

    # ---- ELEMENT GROUPS ----
    REF_ELEM = {}

    tags = KCONEC[9, :]   # 10th row in MATLAB

    REF_ELEM["IN"]          = np.where(tags == 1)[0]
    REF_ELEM["OUT"]         = np.where(tags == 2)[0]
    REF_ELEM["FS"]          = np.where(tags == 4)[0]
    REF_ELEM["BOTTOM"]      = np.where(tags == 5)[0]
    REF_ELEM["IMPERMEABLE"] = np.where(tags == 6)[0]
    REF_ELEM["FRONT"]       = np.where(tags == 7)[0]
    REF_ELEM["BACK"]        = np.where(tags == 8)[0]

    REF_ELEM["POROUS"]      = np.where(tags == 3)[0]
    
    # POROUS_31 ... POROUS_50
    for val in range(31, 51):
        REF_ELEM[f"POROUS_{val}"] = np.where(tags == val)[0]

    return REF_ELEM, NORMAL
