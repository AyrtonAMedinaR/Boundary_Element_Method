import numpy as np

def MATCHING_NODES(NCONEC, KCONEC_1, KCONEC_2, POS_1, POS_2, NODES_MATCH_R1R2, tol=1e-6):

    STORE_NODE_R1 = []
    STORE_NODE_R2 = []

    tol2 = tol**2  # compare squared distance (faster)

    for row in NODES_MATCH_R1R2:
        elem_R1 = int(row[2])
        elem_R2 = int(row[3])

        nodes_R1 = KCONEC_1[:, elem_R1]
        nodes_R2 = KCONEC_2[:, elem_R2]

        coords_R1 = POS_1[nodes_R1]
        coords_R2 = POS_2[nodes_R2]

        # distance matrix between nodes of the two elements
        # shape: (NCONEC, NCONEC)
        diff = coords_R1[:, None, :] - coords_R2[None, :, :]
        dist2 = np.sum(diff**2, axis=2)

        # greedy bijective matching
        used_R2 = set()

        for k1 in range(NCONEC):
            # find closest node in R2
            k2 = np.argmin(dist2[k1])

            if dist2[k1, k2] > tol2:
                raise ValueError(
                    f"Nodes too far apart between elements {elem_R1} and {elem_R2}. "
                    f"Min distance = {np.sqrt(dist2[k1,k2])}"
                )

            if k2 in used_R2:
                raise ValueError(
                    f"Duplicate node match in elements {elem_R1}-{elem_R2}"
                )

            used_R2.add(k2)

            STORE_NODE_R1.append([nodes_R1[k1], elem_R1, k1])
            STORE_NODE_R2.append([nodes_R2[k2], elem_R2, k2])

    STORE_NODE_R1 = np.array(STORE_NODE_R1, dtype=int)
    STORE_NODE_R2 = np.array(STORE_NODE_R2, dtype=int)

    expected = len(NODES_MATCH_R1R2) * NCONEC
    if len(STORE_NODE_R1) != expected:
        raise ValueError(f"Expected {expected} matches, got {len(STORE_NODE_R1)}")

    return STORE_NODE_R1, STORE_NODE_R2