def LINEAR_BC_COMPLEX_G(A, G, NCONEC,
                             STORE_NODE_R1_POS, STORE_NODE_R5_CYL,
                             N0, NE0, N04, NE04, Gporous, k):

    coef = 1j * (Gporous * k)

    # Offsets
    NODE_R1 = STORE_NODE_R1_POS[:, 0] + N0
    ELEM_R1 = STORE_NODE_R1_POS[:, 1] + NE0
    POSI_R1 = STORE_NODE_R1_POS[:, 2]

    NODE_R2 = STORE_NODE_R5_CYL[:, 0] + N04
    ELEM_R2 = STORE_NODE_R5_CYL[:, 1] + NE04
    POSI_R2 = STORE_NODE_R5_CYL[:, 2]

    cols1 = ELEM_R1 * NCONEC + POSI_R1
    cols2 = ELEM_R2 * NCONEC + POSI_R2

    diff12 = G[:, cols1] - G[:, cols2]
    diff21 = -diff12

    A[:, NODE_R1] -= coef * diff12
    A[:, NODE_R2] -= coef * diff21

    return A
