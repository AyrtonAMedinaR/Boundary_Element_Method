import numpy as np

def BC_SCATTERING_OBLIQUE(BC_SCAT, Wave_height, k, omega, gravity, h, THETA, POS, NORMAL, NCONEC, ELEM_SCATT, KCONEC, N, NE, Region):

    # Wave_height: Wave height
    # k: Wavenumber 
    # omega: Frequency
    # gravity: Gravitational accelaration
    # h: Water depth
    # Theta: Angle in radians    
    # POS: x, y and z location of each node
    # NCONEC: Number of nodes in a quad element 
    # ELEM_SCATT: Elements belonging to the region where incident velocity potential is applied
    # KCONEC: Element connectivity array    
    # N: Total number of nodes
    # NE: Total number of elements
    # Region: Region where the BC is applied
    
    N_prev  = sum(N[0:Region]) 
    NE_prev = sum(NE[0:Region])

    kx = np.cos(THETA)
    ky = np.sin(THETA)    

    factor = (gravity * Wave_height / omega)

    for J in ELEM_SCATT + NE_prev:

        nodes = KCONEC[:NCONEC, J]
        XM = POS[nodes - N_prev, 0]
        YM = POS[nodes - N_prev, 1]
        ZM = POS[nodes - N_prev, 2]

        kn = k * ( kx * NORMAL[J - NE_prev,0] + ky * NORMAL[J - NE_prev,1] )

        BC_SCAT[:, J] = (
            factor * kn
            * np.exp(1j * k * ( XM * np.cos(THETA) + YM * np.sin(THETA) ) )
            * (np.cosh(k * (ZM + h)) / np.cosh(k * h))
        )

    return BC_SCAT
