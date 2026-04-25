import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

def VISUAL_FREE_SURFACE(Coef, k, omega, g, Wave_height, hd, NCONEC, ELEM_FS, POS_0, POS_1, KCONEC, N, NE, Region_0):

    from .NODES_FS import NODES_FS

    N_prev = sum(N[0:Region_0])
    
    # NODES AT THE FREE SURFACE
    NODOS_SI = NODES_FS(NCONEC, ELEM_FS, KCONEC, N, NE, Region_0)

    NODOS_UNI = np.unique(NODOS_SI)
    
    PHI = Coef[NODOS_UNI + N_prev]

    xx_BC = POS_0[NODOS_UNI, 0]
    yy_BC = POS_0[NODOS_UNI, 1]   

    zz_BC = np.real((1j * omega / g) * PHI * np.exp(-1j * omega * np.pi / omega))

    # GRID
    xv = np.linspace(np.min(xx_BC), np.max(xx_BC), 501)
    yv = np.linspace(np.min(yy_BC), np.max(yy_BC), 501)

    X_GRID, Y_GRID = np.meshgrid(xv, yv)

    Z_GRID = griddata((xx_BC, yy_BC), zz_BC, (X_GRID, Y_GRID), method='linear')
    
    # Optional scaling for visualization
    scale_z = 1
    Z_plot = Z_GRID * scale_z
    
    # contourLevels = np.linspace(np.nanmin(Z_plot), np.nanmax(Z_plot), 101)
    
    # DOMAIN LIMITS
    x_min, x_max = np.min(POS_0[:,0]), np.max(POS_0[:,0])
    y_min, y_max = np.min(POS_0[:,1]), np.max(POS_0[:,1])
    z_min, z_max = np.min(POS_0[:,2]), np.max(POS_0[:,2])
    
    Lx = x_max - x_min
    Ly = y_max - y_min
    Lz = z_max - z_min
    
    if Lz == 0:
        Lz = 1e-6
    
    # ============================================================
    # FIGURE 1: 3D VIEW
    # ============================================================
    fig1 = plt.figure(figsize=(10,6))
    ax = fig1.add_subplot(111, projection='3d')
    
    surf = ax.plot_surface(X_GRID, Y_GRID, Z_plot, cmap="jet", alpha=0.9, 
                       linewidth=0, 
                       antialiased=False)
    
    ax.scatter(POS_1[:,0], POS_1[:,1], POS_1[:,2],
                color=(0.10,0.10,0.10), s=0.01)
    
    ax.contour(X_GRID, Y_GRID, Z_plot) # , levels=contourLevels)
    
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    ax.set_title("3D Free Surface")
    
    # Aspect ratio
    ax.set_box_aspect([Lx, Ly, Lz])
    
    ax.set_xlim(x_min, x_max)
    ax.set_ylim(y_min, y_max)
    ax.set_zlim(z_min, z_max)
    
    cb1 = fig1.colorbar(surf, ax=ax)
    cb1.set_label(r'$\eta$')
    
    fig1.savefig("free_surface_3D.png", dpi=300, bbox_inches="tight")
    
    # ------------------------------------------------------------
    plt.show()

    # ============================================================
    # FIGURE 2: CENTERLINE
    # ============================================================    
    # index of centerline in y
    j_mid = Y_GRID.shape[0] // 2
    
    x_line = X_GRID[j_mid, :]
    z_line = Z_plot[j_mid, :]
    
    fig1, ax = plt.subplots()
    
    # free surface line
    ax.plot(x_line, z_line, 'r-')
    
    # project structure nodes onto the same plane (optional)
    ax.scatter(POS_1[:,0], POS_1[:,2], s=1)
    
    ax.set_xlabel("x")
    ax.set_ylabel(r'$\eta$')
    ax.set_title("Free Surface (y = mid-plane)")
    
    ax.set_xlim(x_min, x_max)
    
    ax.grid(True)

    #plt.axis('equal')
    
    fig1.savefig("free_surface_2D.png", dpi=300, bbox_inches="tight")
    
    plt.show()

    return zz_BC