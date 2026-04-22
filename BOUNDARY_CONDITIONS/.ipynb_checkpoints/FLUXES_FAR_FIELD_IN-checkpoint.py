import numpy as np

def FLUXES_FAR_FIELD_IN(PHI, DPHI, k, NCONEC, NORMAL, ELEM_SCATT, KCONEC, N, NE, Region):

    # PHI: velocity potential
    # k: wavenumber
    # NCONEC: Number of nodes in a quad element 
    # ELEM_SCATT: Elements belonging to the region where incident velocity potential is applied
    # KCONEC: Element connectivity array    
    # N: Total number of nodes in each region
    # NE: Total number of elements in each region
    # Region: Region where the BC is applied
    
    N_prev  = sum(N[0:Region]) 
    NE_prev = sum(NE[0:Region]) 

    for J in ELEM_SCATT:

        nodes = KCONEC[:NCONEC, J] + N_prev

        kn = k * ( NORMAL[J - NE_prev,0] )

        factor = -(1j * kn)

        DPHI[:, J + NE_prev] = ( factor * PHI[nodes] )

    return DPHI
