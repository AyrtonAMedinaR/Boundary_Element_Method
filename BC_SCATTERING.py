import numpy as np

def BC_SCATTERING(Wave_height, k, omega, h,
                       corner, KCONEC, COLUMNS_R0_LFF,
                       NCONEC, NE0, N0, gravity=9.81):

    BC_SCAT = np.zeros_like(KCONEC, dtype=np.complex128)

    factor = -(gravity * k * Wave_height / omega)

    for J in COLUMNS_R0_LFF + NE0:

        nodes = KCONEC[:, J]
        XM = corner[nodes - N0, 0]
        ZM = corner[nodes - N0, 2]

        BC_SCAT[:, J] = (
            factor
            * np.exp(1j * k * XM)
            * (np.cosh(k * (ZM + h)) / np.cosh(k * h))
        )

    return BC_SCAT
