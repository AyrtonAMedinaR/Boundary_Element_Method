import numpy as np

def POROUS_NODES(NCONEC, KCONEC_1, KCONEC_2, POS_1, POS_2, NODES_MATCH):

    # Pre-round once
    c1 = np.round(POS_1, 4)
    c2 = np.round(POS_2, 4)

    store_r1 = []
    store_r2 = []

    # ---------------------------------------------------------
    # Build lookup table for REGION 2
    # key = (x,y,z)  -> [(NODE2, J2, K2), ...]
    # ---------------------------------------------------------
    lookup = {}

    for J2 in NODES_MATCH[:, 3]:
        nodes2 = KCONEC_2[:, J2]

        for K2, NODE2 in enumerate(nodes2):
            key = tuple(c2[NODE2])
            lookup.setdefault(key, []).append((NODE2, J2, K2))

    # ---------------------------------------------------------
    # Search matches for REGION 1
    # ---------------------------------------------------------
    for J1 in NODES_MATCH[:, 2]:
        nodes1 = KCONEC_1[:, J1]

        for K1, NODE1 in enumerate(nodes1):
            key = tuple(c1[NODE1])

            if key in lookup:
                for (NODE2, J2, K2) in lookup[key]:
                    store_r1.append([NODE1, J1, K1])
                    store_r2.append([NODE2, J2, K2])

    STORE_1 = np.asarray(store_r1, dtype=int)
    STORE_2 = np.asarray(store_r2, dtype=int)

    # ---------------------------------------------------------
    # Remove consecutive duplicates (vectorized)
    # ---------------------------------------------------------
    if len(STORE_1) > 0:
        mask = np.ones(len(STORE_1), dtype=bool)
        mask[1:] = np.any(STORE_1[1:] != STORE_1[:-1], axis=1)
        STORE_1 = STORE_1[mask]
        STORE_2 = STORE_2[mask]

    return STORE_1, STORE_2


    