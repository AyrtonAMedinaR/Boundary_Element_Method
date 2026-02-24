import numpy as np

def POROUS_NODES(NCONEC, KCONEC1, KCONEC2, corner1, corner2, NODES_MATCH_R1R2):

    # Pre-round once
    c1 = np.round(corner1, 4)
    c2 = np.round(corner2, 4)

    store_r1 = []
    store_r2 = []

    # ---------------------------------------------------------
    # Build lookup table for REGION 2
    # key = (x,y,z)  -> [(NODE2, J2, K2), ...]
    # ---------------------------------------------------------
    lookup = {}

    for J2 in NODES_MATCH_R1R2[:, 3]:
        nodes2 = KCONEC2[:, J2]

        for K2, NODE2 in enumerate(nodes2):
            key = tuple(c2[NODE2])
            lookup.setdefault(key, []).append((NODE2, J2, K2))

    # ---------------------------------------------------------
    # Search matches for REGION 1
    # ---------------------------------------------------------
    for J1 in NODES_MATCH_R1R2[:, 2]:
        nodes1 = KCONEC1[:, J1]

        for K1, NODE1 in enumerate(nodes1):
            key = tuple(c1[NODE1])

            if key in lookup:
                for (NODE2, J2, K2) in lookup[key]:
                    store_r1.append([NODE1, J1, K1])
                    store_r2.append([NODE2, J2, K2])

    R1 = np.asarray(store_r1, dtype=int)
    R2 = np.asarray(store_r2, dtype=int)

    # ---------------------------------------------------------
    # Remove consecutive duplicates (vectorized)
    # ---------------------------------------------------------
    if len(R1) > 0:
        mask = np.ones(len(R1), dtype=bool)
        mask[1:] = np.any(R1[1:] != R1[:-1], axis=1)
        R1 = R1[mask]
        R2 = R2[mask]

    return R1, R2


    