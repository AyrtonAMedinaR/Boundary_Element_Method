import numpy as np

def FLUXES_FAR_FIELD_IN(PHI, DPHI, k, NCONEC, NORMAL, ELEM_SCATT, KCONEC, N, NE, Region):

    # PHI: array with velocity potentials
    # DPHI: array with fluxes    
    # k: wavenumber
    # NCONEC: Number of nodes in a quad element 
    # NORMAL: Array containing the normal vectors of each element    
    # ELEM_SCATT: Elements belonging to the region where incident velocity potential is applied
    # KCONEC: Element connectivity array    
    # N: Total number of nodes
    # NE: Total number of elements
    # Region: Region where the BC is applied
    
    N_prev  = sum(N[0:Region]) 
    NE_prev = sum(NE[0:Region]) 

    for J in ELEM_SCATT:

        nodes = KCONEC[:NCONEC, J] + N_prev

        kn = k * ( NORMAL[J - NE_prev,0] )

        factor = -(1j * kn)

        DPHI[:, J + NE_prev] = ( factor * PHI[nodes] )

    return DPHI
