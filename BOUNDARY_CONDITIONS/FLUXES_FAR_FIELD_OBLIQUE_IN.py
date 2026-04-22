import numpy as np

def FLUXES_FAR_FIELD_OBLIQUE_IN(PHI, DPHI, k, THETA, NCONEC, NORMAL, ELEM_SCATT, KCONEC, N, NE, Region):

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

    kx = np.cos(THETA)
    ky = np.sin(THETA)    

    # for J in ELEM_SCATT + NE_prev:
    for J in ELEM_SCATT:

        # nodes = KCONEC[:, J]
        # nodes = KCONEC[0:NCONEC, J]
        nodes = KCONEC[:NCONEC, J] + N_prev

        kn = k * ( kx * NORMAL[J - NE_prev,0] + ky * NORMAL[J - NE_prev,1] )

        factor = -(1j * k)

        DPHI[:, J + NE_prev] = ( factor * PHI[nodes] )

    return DPHI
