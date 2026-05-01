# matching_interface.py

# Utilities to match interface elements and nodes between two meshes.

# Workflow:
# 1) MATCHING(...) -> finds matching interface ELEMENTS via center node.
# 2) MATCHING_NODES(...) -> finds matching NODES inside each matched element pair.
# 3) MATCH_INTERFACES(...) -> convenience wrapper that runs both steps.

import numpy as np
from scipy.spatial import cKDTree

# ==========================================================
# STEP 1 — MATCH INTERFACE ELEMENTS (via center node)
# ==========================================================
def MATCHING(REF_1, POS_1, KCONEC_1,
             REF_2, POS_2, KCONEC_2,
             tol=1e-6):
    
    # Match interface elements between region 1 and region 2 using the
    # element center node (node index 8).

    # Returns
    # NODES_MATCH : (N,4) int array
    #     col 0 → center node in R1
    #     col 1 → center node in R2
    #     col 2 → element index in R1
    #     col 3 → element index in R2

    # ---- interface center nodes ----
    nodes_R1 = KCONEC_1[8, REF_1]
    nodes_R2 = KCONEC_2[8, REF_2]

    coords_R1 = POS_1[nodes_R1]
    coords_R2 = POS_2[nodes_R2]

    # ---- KD-tree nearest search ----
    tree = cKDTree(coords_R2)
    dist, idx = tree.query(coords_R1, distance_upper_bound=tol)

    if np.any(np.isinf(dist)):
        bad = np.where(np.isinf(dist))[0][0]
        raise ValueError(
            "No matching interface node found within tolerance. "
            f"Distance = {dist[bad]}"
        )

    # ---- build mapping table ----
    N = len(nodes_R1)
    NODES_MATCH = np.zeros((N, 4), dtype=int)

    NODES_MATCH[:, 0] = nodes_R1
    NODES_MATCH[:, 1] = nodes_R2[idx]

    # fast node → element maps
    node_to_elem_R1 = {KCONEC_1[8, e]: e for e in REF_1}
    node_to_elem_R2 = {KCONEC_2[8, e]: e for e in REF_2}

    NODES_MATCH[:, 2] = [node_to_elem_R1[n] for n in nodes_R1]
    NODES_MATCH[:, 3] = [node_to_elem_R2[n] for n in NODES_MATCH[:, 1]]

    return NODES_MATCH
    
# ==========================================================
# STEP 2 — MATCH NODES INSIDE EACH ELEMENT PAIR
# ==========================================================
def MATCHING_NODES(NCONEC, KCONEC_1, KCONEC_2,
                   POS_1, POS_2,
                   NODES_MATCH_R1R2,
                   tol=1e-6):
    
    # For each matched element pair, match all local nodes.

    # Returns
    # STORE_NODE_R1 : (Npairs*NCONEC,3)
    # STORE_NODE_R2 : (Npairs*NCONEC,3)

    # Each row: [global_node_id, element_id, local_node_index]

    STORE_NODE_R1 = []
    STORE_NODE_R2 = []
    tol2 = tol**2

    for row in NODES_MATCH_R1R2:
        elem_R1 = int(row[2])
        elem_R2 = int(row[3])

        nodes_R1 = KCONEC_1[:, elem_R1]
        nodes_R2 = KCONEC_2[:, elem_R2]

        coords_R1 = POS_1[nodes_R1]
        coords_R2 = POS_2[nodes_R2]

        # pairwise distance matrix
        diff = coords_R1[:, None, :] - coords_R2[None, :, :]
        dist2 = np.sum(diff**2, axis=2)

        used_R2 = set()

        for k1 in range(NCONEC):
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


# ==========================================================
# OPTIONAL — FULL PIPELINE WRAPPER
# ==========================================================
def MATCHING_INTERFACES(REF_1, POS_1, KCONEC_1,
                     REF_2, POS_2, KCONEC_2,
                     NCONEC,
                     tol=1e-5):
    
    # Convenience wrapper that performs the full matching pipeline.

    NODES_MATCH = MATCHING(
        REF_1, POS_1, KCONEC_1,
        REF_2, POS_2, KCONEC_2,
        tol
    )

    STORE_NODE_R1, STORE_NODE_R2 = MATCHING_NODES(
        NCONEC,
        KCONEC_1, KCONEC_2,
        POS_1, POS_2,
        NODES_MATCH,
        tol
    )

    return STORE_NODE_R1, STORE_NODE_R2