import numpy as np
from scipy.interpolate import griddata

def ENERGY_FLUX(Coef, NCONEC, omega, rho, k, g, Wave_height, hd,
                    ELEM_R, POS_R, KCONEC_R,
                    ELEM_T, POS_T, KCONEC_T):

    from NODES_FS import NODES_FS

    # -------------------------------------------------------
    # NODES AT FREE SURFACE
    # -------------------------------------------------------
    NODOS_R = NODES_FS(NCONEC, KCONEC_R, ELEM_R, 0)
    NODOS_T = NODES_FS(NCONEC, KCONEC_T, ELEM_T, 0)

    NODOS_UNI_R = np.unique(NODOS_R)
    NODOS_UNI_T = np.unique(NODOS_T)

    # -------------------------------------------------------
    # POTENTIALS
    # -------------------------------------------------------
    PHI_R = (
        Coef[NODOS_UNI_R]
        - (-1j * g * Wave_height / (2 * omega))
        * (np.cosh(k * (POS_R[NODOS_UNI_R, 2] + hd)) / np.cosh(k * hd))
        * np.exp(1j * POS_R[NODOS_UNI_R, 0])
    )

    PHI_T = Coef[NODOS_UNI_T]

    # -------------------------------------------------------
    # ENERGY FLUX DENSITIES
    # -------------------------------------------------------
    FLUX_R = k * omega * rho * np.abs(PHI_R)**2
    FLUX_T = k * omega * rho * np.abs(PHI_T)**2

    COS_VARTHETA_R = np.cos(0.0)
    COS_VARTHETA_T = -np.cos(np.pi)

    # -------------------------------------------------------
    # COORDINATES
    # -------------------------------------------------------
    xx_R = POS_R[NODOS_UNI_R, 1]
    yy_R = POS_R[NODOS_UNI_R, 2]
    zz_R = (FLUX_R * COS_VARTHETA_R) / 2

    xx_T = POS_T[NODOS_UNI_T, 1]
    yy_T = POS_T[NODOS_UNI_T, 2]
    zz_T = (FLUX_T * COS_VARTHETA_T) / 2

    # -------------------------------------------------------
    # GRID INTERPOLATION
    # -------------------------------------------------------
    xv_R = np.linspace(xx_R.min(), xx_R.max(), 101)
    yv_R = np.linspace(yy_R.min(), yy_R.max(), 101)
    X_R, Y_R = np.meshgrid(xv_R, yv_R)
    Z_R = griddata((xx_R, yy_R), zz_R, (X_R, Y_R), method="linear")

    xv_T = np.linspace(xx_T.min(), xx_T.max(), 101)
    yv_T = np.linspace(yy_T.min(), yy_T.max(), 101)
    X_T, Y_T = np.meshgrid(xv_T, yv_T)
    Z_T = griddata((xx_T, yy_T), zz_T, (X_T, Y_T), method="linear")

    # Replace NaNs produced by griddata
    Z_R = np.nan_to_num(Z_R)
    Z_T = np.nan_to_num(Z_T)

    # -------------------------------------------------------
    # NUMERICAL INTEGRATION
    # -------------------------------------------------------
    FLUX_REFLE = np.trapezoid(np.trapezoid(Z_R, xv_R, axis=1), yv_R)
    FLUX_TRANS = np.trapezoid(np.trapezoid(Z_T, xv_T, axis=1), yv_T)

    return FLUX_REFLE, FLUX_TRANS
