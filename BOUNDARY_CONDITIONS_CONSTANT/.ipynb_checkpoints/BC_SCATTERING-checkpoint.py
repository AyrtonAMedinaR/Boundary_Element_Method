import numpy as np

def BC_SCATTERING(BC_SCAT, Wave_height, k, omega, gravity, h, POS, ELEM_SCATT, N, NE, Region):

    # Wave_height: Wave height
    # k: Wavenumber 
    # omega: Frequency
    # gravity: Gravitational accelaration
    # h: Water depth
    # POS: x, y and z location of each node
    # NCONEC: Number of nodes in a quad element 
    # ELEM_SCATT: Elements belonging to the region where incident velocity potential is applied
    # KCONEC: Element connectivity array    
    # N: Total number of nodes in each region
    # NE: Total number of elements in each region
    # Region: Region where the BC is applied
    
    N_prev  = sum(N[0:Region]) 
    NE_prev = sum(NE[0:Region]) 

    factor = -(gravity * k * Wave_height / omega)

    for J in ELEM_SCATT + NE_prev:

        XM = POS[J - N_prev, 0]
        ZM = POS[J - N_prev, 2]

        BC_SCAT[J] = (
            factor
            * np.exp(1j * k * XM)
            * (np.cosh(k * (ZM + h)) / np.cosh(k * h))
        )

    return BC_SCAT
