import numpy as np

def BC_FREE_SURFACE(A, B, omega, g, ELEM_FS, NE, Region):

    # A and B: Coefficient matrices 
    # omega: Frequency 
    # gravity: Gravitational acceleration 
    # NCONEC: Number of nodes in a quad element 
    # ELEM_FS: Elements belonging to the free surface
    # KCONEC: Element connectivity array
    # NE: Total number of elements in each region
    # Region: Region where the BC is applied

    NE_prev = sum(NE[0:Region])
    
    factor = (omega**2) / g

    for J in ELEM_FS + NE_prev:
        
            A[:, J] = A[:, J] - factor * B[:, J]
  
    return A
