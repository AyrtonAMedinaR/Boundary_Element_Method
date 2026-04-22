import numpy as np

def BC_FAR_FIELD_OUT(A, B, k, NCONEC, NORMAL, ELEM_FF, KCONEC, NE, Region):

    # A and B: Coefficient matrices 
    # k: wavenumber 
    # NCONEC: Number of nodes in a quad element 
    # NORMAL: Array [NEx3] containing the normal vectors of each element
    # ELEM_FF: Elements belonging to the far field
    # KCONEC: Element connectivity array
    # NE: Total number of elements in each region
    # Region: Region where the BC is applied  
    
    NE_prev = sum(NE[0:Region]) 

    for J in ELEM_FF + NE_prev:
     
        IKs = KCONEC[:NCONEC, J]     
        cols_B = J*NCONEC + np.arange(NCONEC)

        kn = k * ( NORMAL[J - NE_prev,0] )

        A[:, IKs] -= 1j * kn * B[:, cols_B]     

    return A

