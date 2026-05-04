import numpy as np
from scipy.interpolate import griddata

def INT_FLUX_FS_S(Coef, omega, g, NCONEC, ELEM, POS, KCONEC, N, NE, Region):

    # Coef: Vector containing PHI values  
    # omega: Frequency 
    # gravity: Gravitational acceleration     
    # NCONEC: Number of nodes in a quad element 
    # ELEM: Elements belonging to this BC
    # POS: x, y and z location of each node    
    # KCONEC: Element connectivity array    
    # N: Total number of nodes
    # NE: Total number of elements
    # Region: Region where the BC is applied
    
    from BEM_UTILS.NODES_FS import NODES_FS

    N_prev  = sum(N[0:Region])  

    # NODES AT BOUNDARY
    NODOS = NODES_FS(NCONEC, ELEM, KCONEC, N, NE, Region)

    NODOS_UNI = np.unique(NODOS) + N_prev

    PHI_R = Coef[NODOS_UNI]

    xx_R1 = POS[NODOS_UNI - N_prev, 0]
    yy_R1 = POS[NODOS_UNI - N_prev, 1]
    zz_R1 = (omega**2 / g) * PHI_R

    # --- GRID ---
    xv_R1 = np.linspace(np.min(xx_R1), np.max(xx_R1), 101)
    yv_R1 = np.linspace(np.min(yy_R1), np.max(yy_R1), 101)

    X_R1, Y_R1 = np.meshgrid(xv_R1, yv_R1)

    # Interpolation (MATLAB griddata equivalent)
    Z_R1 = griddata( (xx_R1, yy_R1), zz_R1, (X_R1, Y_R1), method='linear' )

    # Replace NaN with 0 (circular domain → square grid)
    Z_R1 = np.nan_to_num(Z_R1, nan=0.0)

    # --- DOUBLE INTEGRAL using trapezoidal rule ---
    FLUX = np.trapezoid( np.trapezoid(Z_R1, xv_R1, axis=1), yv_R1, axis=0 )

    return FLUX