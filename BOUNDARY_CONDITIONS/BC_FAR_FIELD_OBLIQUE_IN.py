import numpy as np

def BC_FAR_FIELD_OBLIQUE_IN(A, B, k, THETA, NCONEC, NORMAL, ELEM_FF, KCONEC, NE, Region):

    # A and B: Coefficient matrices 
    # k: Wavenumber 
    # Theta: Angle in radians
    # NCONEC: Number of nodes in a quad element 
    # NORMAL: Array [NEx3] containing the normal vectors of each element
    # ELEM_FF: Elements belonging to the far field
    # KCONEC: Element connectivity array
    # NE: Total number of elements
    # Region: Region where the BC is applied  

    kx = np.cos(THETA)
    ky = np.sin(THETA)
    
    NE_prev = sum(NE[0:Region]) 

    for J in ELEM_FF + NE_prev:

        # IKs = KCONEC[:, J]     
        IKs = KCONEC[0:NCONEC, J]     
        cols_B = J*NCONEC + np.arange(NCONEC)

        kn = k * ( kx * NORMAL[J - NE_prev,0] + ky * NORMAL[J - NE_prev,1] )

        A[:, IKs] += 1j * kn * B[:, cols_B]

    return A

