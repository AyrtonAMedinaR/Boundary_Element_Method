import numpy as np

def FLUXES_FREE_SURFACE(PHI, DPHI, omega, gravity, NCONEC, ELEM_SCATT, KCONEC, N, NE, Region):

    # PHI: velocity potential
    # omega: Frequency
    # gravity: Gravitational accelaration
    # NCONEC: Number of nodes in a quad element 
    # ELEM_SCATT: Elements belonging to the region where incident velocity potential is applied
    # KCONEC: Element connectivity array    
    # N: Total number of nodes in each region
    # NE: Total number of elements in each region
    # Region: Region where the BC is applied
    
    N_prev  = sum(N[0:Region]) 
    NE_prev = sum(NE[0:Region]) 

    factor = (omega**2/gravity)

    for J in ELEM_SCATT + NE_prev:

        nodes = KCONEC[:, J]

        DPHI[:, J] = ( factor * PHI[nodes - N_prev] )

    return DPHI
