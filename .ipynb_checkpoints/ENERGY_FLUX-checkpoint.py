import numpy as np
from scipy.interpolate import griddata

def ENERGY_FLUX(Coef, NCONEC, omega, rho, k, g, Wave_height, hd,
                    COLUMNS_R0_LFF, corner0, KCONEC0,
                    COLUMNS_R6_RFF, corner6, KCONEC6,
                    N05):

    from NODES_FS import NODES_FS

    # -------------------------------------------------------
    # NODES AT FREE SURFACE
    # -------------------------------------------------------
    NODOS_R0_SI = NODES_FS(NCONEC, KCONEC0, COLUMNS_R0_LFF, 0)
    NODOS_R6_SI = NODES_FS(NCONEC, KCONEC6, COLUMNS_R6_RFF, 0)

    NODOS_SF_R0 = np.unique(NODOS_R0_SI)
    NODOS_SF_R6 = np.unique(NODOS_R6_SI) + N05

    # -------------------------------------------------------
    # POTENTIALS
    # -------------------------------------------------------
    PHI_R = (
        Coef[NODOS_SF_R0]
        - (-1j * g * Wave_height / (2 * omega))
        * (np.cosh(k * (corner0[NODOS_SF_R0, 2] + hd)) / np.cosh(k * hd))
        * np.exp(1j * corner0[NODOS_SF_R0, 0])
    )

    PHI_T = Coef[NODOS_SF_R6]

    # -------------------------------------------------------
    # ENERGY FLUX DENSITIES
    # -------------------------------------------------------
    FLUX_R0 = k * omega * rho * np.abs(PHI_R)**2
    FLUX_R6 = k * omega * rho * np.abs(PHI_T)**2

    COS_VARTHETA_R = np.cos(0.0)
    COS_VARTHETA_T = -np.cos(np.pi)

    # -------------------------------------------------------
    # COORDINATES
    # -------------------------------------------------------
    xx_R0 = corner0[NODOS_SF_R0, 1]
    yy_R0 = corner0[NODOS_SF_R0, 2]
    zz_R0 = (FLUX_R0 * COS_VARTHETA_R) / 2

    idx6 = NODOS_SF_R6 - N05
    xx_R6 = corner6[idx6, 1]
    yy_R6 = corner6[idx6, 2]
    zz_R6 = (FLUX_R6 * COS_VARTHETA_T) / 2

    # -------------------------------------------------------
    # GRID INTERPOLATION
    # -------------------------------------------------------
    xv_R0 = np.linspace(xx_R0.min(), xx_R0.max(), 101)
    yv_R0 = np.linspace(yy_R0.min(), yy_R0.max(), 101)
    X_R0, Y_R0 = np.meshgrid(xv_R0, yv_R0)
    Z_R0 = griddata((xx_R0, yy_R0), zz_R0, (X_R0, Y_R0), method="linear")

    xv_R6 = np.linspace(xx_R6.min(), xx_R6.max(), 101)
    yv_R6 = np.linspace(yy_R6.min(), yy_R6.max(), 101)
    X_R6, Y_R6 = np.meshgrid(xv_R6, yv_R6)
    Z_R6 = griddata((xx_R6, yy_R6), zz_R6, (X_R6, Y_R6), method="linear")

    # Replace NaNs produced by griddata
    Z_R0 = np.nan_to_num(Z_R0)
    Z_R6 = np.nan_to_num(Z_R6)

    # -------------------------------------------------------
    # NUMERICAL INTEGRATION
    # -------------------------------------------------------
    FLUX_REFLE = np.trapezoid(np.trapezoid(Z_R0, xv_R0, axis=1), yv_R0)
    FLUX_TRANS = np.trapezoid(np.trapezoid(Z_R6, xv_R6, axis=1), yv_R6)

    return FLUX_REFLE, FLUX_TRANS
