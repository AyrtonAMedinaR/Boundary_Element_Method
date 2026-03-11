import numpy as np

def BC_SCATTERING(Wave_height, k, omega, gravity, h, corner, KCONEC, COLUMNS_R0_LFF, NCONEC, N, NE, Region):

    N_prev  = sum(N[0:Region]) 
    NE_prev = sum(NE[0:Region]) 

    BC_SCAT = np.zeros_like(KCONEC, dtype=np.complex128)

    factor = -(gravity * k * Wave_height / omega)

    for J in COLUMNS_R0_LFF + NE_prev:

        nodes = KCONEC[:, J]
        XM = corner[nodes - N_prev, 0]
        ZM = corner[nodes - N_prev, 2]

        BC_SCAT[:, J] = (
            factor
            * np.exp(1j * k * XM)
            * (np.cosh(k * (ZM + h)) / np.cosh(k * h))
        )

    return BC_SCAT
