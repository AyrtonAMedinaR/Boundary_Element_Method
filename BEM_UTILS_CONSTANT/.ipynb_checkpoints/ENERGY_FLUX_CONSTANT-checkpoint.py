import numpy as np
from scipy.interpolate import griddata

def ENERGY_FLUX_CONSTANT(Coef, omega, rho, k, g, Wave_height, hd, AREA_ELEM,
                ELEM_R, POS_R,
                ELEM_T, POS_T,
                N, NE, Region_R, Region_T):

    # Indexing logic for the specific regions
    N_prev_R = sum(N[0:Region_R]) 
    N_prev_T = sum(N[0:Region_T])  

    # Mapping to global coefficient indices
    NODOS_UNI_R = ELEM_R + N_prev_R
    NODOS_UNI_T = ELEM_T + N_prev_T

    # POTENTIALS at the centroid position
    PHI_R = (
        Coef[NODOS_UNI_R]
        - (-1j * g * Wave_height / (2 * omega))
        * (np.cosh(k * (POS_R[ELEM_R, 2] + hd)) / np.cosh(k * hd))
        * np.exp(1j * k * POS_R[ELEM_R, 0]) # Added 'k' here for spatial phase
    )

    PHI_T = Coef[NODOS_UNI_T]

    # ENERGY FLUX DENSITIES
    FLUX_R = k * omega * rho * np.abs(PHI_R)**2
    FLUX_T = k * omega * rho * np.abs(PHI_T)**2

    # Normalization factors
    COS_VARTHETA_R = np.cos(0.0)
    COS_VARTHETA_T =-np.cos(np.pi)

    # NUMERICAL INTEGRATION (Summation over Area)
    # We multiply the flux density by the area of each constant element
    FLUX_REFLE = np.sum(FLUX_R * AREA_ELEM[NODOS_UNI_R]) * COS_VARTHETA_R / 2
    FLUX_TRANS = np.sum(FLUX_T * AREA_ELEM[NODOS_UNI_T]) * COS_VARTHETA_T / 2

    return FLUX_REFLE, FLUX_TRANS