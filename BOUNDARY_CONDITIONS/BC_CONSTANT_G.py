def BC_CONSTANT_G(A, B, DEPTH, NCONEC, R0_STRUCT, R1_STRUCT, N, NE, Gporous, Region_M, Region_N):

    # A and B: Coefficient matrices 
    # k: Wavenumber 
    # NCONEC: Number of nodes in a quad element 
    # R0_STRUCT: Nodes, elements and position (from 1 to 9) of structure in region 0
    # R1_STRUCT: Nodes, elements and position (from 1 to 9) of structure in region 1 
    # N: Total number of nodes in each region
    # NE: Total number of elements in each region
    # Gporous: Porosity parameter
    # Region_M: Outer region where the BC is applied   
    # Region_N: Inner region where the BC is applied   
    
    N_prev_M  = sum(N[0:Region_M])
    N_prev_N  = sum(N[0:Region_N])

    NE_prev_M = sum(NE[0:Region_M])
    NE_prev_N = sum(NE[0:Region_N])
    
    coef = 1j * (Gporous / DEPTH)

    # Offsets
    NODE_R_M = R0_STRUCT[:, 0] + N_prev_M
    ELEM_R_M = R0_STRUCT[:, 1] + NE_prev_M
    POSI_R_M = R0_STRUCT[:, 2]

    NODE_R_N = R1_STRUCT[:, 0] + N_prev_N
    ELEM_R_N = R1_STRUCT[:, 1] + NE_prev_N
    POSI_R_N = R1_STRUCT[:, 2]

    cols1 = ELEM_R_M * NCONEC + POSI_R_M
    cols2 = ELEM_R_N * NCONEC + POSI_R_N

    diff12 = B[:, cols1] - B[:, cols2]
    diff21 = -diff12

    A[:, NODE_R_M] -= coef * diff12
    A[:, NODE_R_N] -= coef * diff21

    return A
