from scipy.spatial import cKDTree
import numpy as np

def MATCHING(REF_1, POS_1, KCONEC_1, REF_2, POS_2, KCONEC_2, tol=1e-6):

    # ---- interface nodes (center nodes of elements) ----
    nodes_R1 = KCONEC_1[8, REF_1]
    nodes_R2 = KCONEC_2[8, REF_2]

    coords_R1 = POS_1[nodes_R1]
    coords_R2 = POS_2[nodes_R2]

    # ---- KD-tree search ----
    tree = cKDTree(coords_R2)
    dist, idx = tree.query(coords_R1, distance_upper_bound=tol)

    # check failures
    if np.any(np.isinf(dist)):
        bad = np.where(np.isinf(dist))[0][0]
        raise ValueError(
            f"No matching interface node found within tol. "
            f"Distance = {dist[bad]}"
        )

    # ---- build NODES_MATCH ----
    N = len(nodes_R1)
    NODES_MATCH = np.zeros((N,4), dtype=int)

    NODES_MATCH[:,0] = nodes_R1                # node in R1
    NODES_MATCH[:,1] = nodes_R2[idx]          # matching node in R2

    # ---- build fast node→element maps ----
    node_to_elem_R1 = {KCONEC_1[8,e]: e for e in REF_1}
    node_to_elem_R2 = {KCONEC_2[8,e]: e for e in REF_2}

    NODES_MATCH[:,2] = [node_to_elem_R1[n] for n in nodes_R1]
    NODES_MATCH[:,3] = [node_to_elem_R2[n] for n in NODES_MATCH[:,1]]

    return NODES_MATCH