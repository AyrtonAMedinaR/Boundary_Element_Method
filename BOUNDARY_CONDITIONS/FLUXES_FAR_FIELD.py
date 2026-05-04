import numpy as np

def FLUXES_FAR_FIELD(PHI, DPHI, k, NCONEC, ELEM_SCATT, KCONEC, N, NE, Region):

    # PHI: velocity potential
    # DPHI: array with fluxes        
    # k: wavenumber
    # NCONEC: Number of nodes in a quad element 
    # ELEM_SCATT: Elements belonging to the region where incident velocity potential is applied
    # KCONEC: Element connectivity array    
    # N: Total number of nodes
    # NE: Total number of elements
    # Region: Region where the BC is applied
    
    N_prev  = sum(N[0:Region]) 
    NE_prev = sum(NE[0:Region]) 

    factor = (1j * k)

    # for J in ELEM_SCATT + NE_prev:
    for J in ELEM_SCATT:

        # nodes = KCONEC[:, J]
        # nodes = KCONEC[0:NCONEC, J]
        nodes = KCONEC[0:NCONEC, J] + N_prev

        DPHI[:, J + NE_prev] = ( factor * PHI[nodes] )

    return DPHI
