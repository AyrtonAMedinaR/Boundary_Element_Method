import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def VISUAL_FREE_SURFACE(Coef, k, omega, g, Wave_height, hd, NCONEC, ELEM_FS, POS_0, POS_1, KCONEC, N, NE, Region_0):

    from .NODES_FS import NODES_FS

    # NODES AT THE FREE SURFACE
    NODOS_SI = NODES_FS(NCONEC, ELEM_FS, KCONEC, N, NE, Region_0)

    NODOS_UNI = np.unique(NODOS_SI)

    PHI = Coef[NODOS_UNI]

    xx_BC = POS_0[NODOS_UNI, 0]
    yy_BC = POS_0[NODOS_UNI, 1]

    zz_BC = np.real((1j * omega / g) * PHI * np.exp(-1j * omega * np.pi / omega)) / (Wave_height / 2)

    # GRID
    xv = np.linspace(np.min(xx_BC), np.max(xx_BC), 501)
    yv = np.linspace(np.min(yy_BC), np.max(yy_BC), 501)

    X_GRID, Y_GRID = np.meshgrid(xv, yv)

    Z_GRID = griddata((xx_BC, yy_BC), zz_BC, (X_GRID, Y_GRID), method='linear')

    # PLOT
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    
    surf = ax.plot_surface(X_GRID, Y_GRID, Z_GRID, cmap="jet", alpha=0.25)
    
    ax.scatter(POS_1[:,0], POS_1[:,1], POS_1[:,2],
               color=(0.05,0.05,0.05), s=2)
    
    contourLevels = np.linspace(-1.5, 1.5, 101)
    
    ax.contour(X_GRID, Y_GRID, Z_GRID, levels=contourLevels)
    
    cb = fig.colorbar(surf, ax=ax)
    cb.set_label(r'$\bar{\eta}$')
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")

    ax.set_box_aspect([7.5,5,1])  # axis equal
    
    # surf.set_clim(-1.5, 1.5)

    fig.savefig("free_surface.png", dpi=300, bbox_inches="tight")
    
    plt.show()

    return PHI