import numpy as np

def ELEMENT_CONNECTION_BREEF(NCONEC, KCONEC1, POS1, NE1):
    # """
    # NCONEC  : number of nodes per element (9)
    # KCONEC1 : (10, NE1) array
    # POS1    : (N1,3) array
    # NE1     : number of elements

    # Returns
    # -------
    # REF_ELEM : dictionary of element groups
    # NORMAL1  : (NE1,3) array of normals
    # """

    from COS_DIR import COS_DIR

    NORMAL1 = np.zeros((NE1, 3))

    # Coordinates
    X = POS1[:,0]
    Y = POS1[:,1]
    Z = POS1[:,2]

    # ---- NORMAL VECTORS ----
    for J in range(NE1):
        ETA1 = COS_DIR(KCONEC1, J, X, Y, Z)

        NORMAL1[J,0] = ETA1[0]
        NORMAL1[J,1] = ETA1[1]
        NORMAL1[J,2] = ETA1[2]

    # ---- ELEMENT GROUPS ----
    REF_ELEM = {}

    tags = KCONEC1[9, :]   # 10th row in MATLAB

    REF_ELEM["IN"]          = np.where(tags == 1)[0]
    REF_ELEM["OUT"]         = np.where(tags == 2)[0]
    REF_ELEM["FS"]          = np.where(tags == 4)[0]
    REF_ELEM["BOTTOM"]      = np.where(tags == 5)[0]
    REF_ELEM["IMPERMEABLE"] = np.where(tags == 6)[0]
    REF_ELEM["FRONT"]       = np.where(tags == 7)[0]
    REF_ELEM["BACK"]        = np.where(tags == 8)[0]

    # POROUS_31 ... POROUS_50
    for val in range(31, 51):
        REF_ELEM[f"POROUS_{val}"] = np.where(tags == val)[0]

    return REF_ELEM, NORMAL1
