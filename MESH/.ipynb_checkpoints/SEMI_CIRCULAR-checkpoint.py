import sys
import gmsh
import numpy as np
import os
import natsort
import math
from pathlib import Path

def setRecombine(SV_params):
    import gmsh
    for localParam in SV_params:
        gmsh.model.mesh.setRecombine(2,localParam)

def SEMI_CIRCULAR(hd,RADIUS,WIDTH,LX):
    
    # NUMBERS OF ELEMENTS FOR MESH
    ELEM = 7;
    ELEM_ARCS = 6;
      
    DX = 0.10*hd;

    # LENGTH OF THE ELEMENTS IN GMSH
    lc = hd/10;
       
    path = os.getcwd()   
    os.chdir(path)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    gmsh.model.add("modelo_semi_circular")

    ###############################
    # The following line allows you to abbreviate the code:
    geom = gmsh.model.geo
    ###############################

    ###############################
    # POINTS
    p = [None]*26
    c = math.cos(math.radians(45))
    
    # -------------------------
    # BOTTOM
    # -------------------------
    p[0] = geom.addPoint(LX, -WIDTH/2, -hd, lc)
    p[1] = geom.addPoint(LX + 2*DX + 2*RADIUS, -WIDTH/2, -hd, lc)
    p[2] = geom.addPoint(LX, WIDTH/2, -hd, lc)
    p[3] = geom.addPoint(LX + 2*DX + 2*RADIUS, WIDTH/2, -hd, lc)
    
    # -------------------------
    # FREE SURFACE
    # -------------------------
    p[4] = geom.addPoint(LX, -WIDTH/2, 0, lc)
    p[5] = geom.addPoint(LX + 2*DX + 2*RADIUS, -WIDTH/2, 0, lc)
    p[6] = geom.addPoint(LX, WIDTH/2, 0, lc)
    p[7] = geom.addPoint(LX + 2*DX + 2*RADIUS, WIDTH/2, 0, lc)
    
    # -------------------------
    # HALF CYLINDER - FRONT
    # -------------------------
    p[8]  = geom.addPoint(LX + DX, -WIDTH/2, -hd, lc)
    p[9]  = geom.addPoint(LX + DX + RADIUS, -WIDTH/2, -hd, lc)
    p[10] = geom.addPoint(LX + DX + 2*RADIUS, -WIDTH/2, -hd, lc)
    p[11] = geom.addPoint(LX + DX + RADIUS - RADIUS*c, -WIDTH/2, -hd + RADIUS*c, lc)
    p[12] = geom.addPoint(LX + DX + RADIUS + RADIUS*c, -WIDTH/2, -hd + RADIUS*c, lc)
    
    # -------------------------
    # HALF CYLINDER - BACK
    # -------------------------
    p[13] = geom.addPoint(LX + DX, WIDTH/2, -hd, lc)
    p[14] = geom.addPoint(LX + DX + RADIUS, WIDTH/2, -hd, lc)
    p[15] = geom.addPoint(LX + DX + 2*RADIUS, WIDTH/2, -hd, lc)
    p[16] = geom.addPoint(LX + DX + RADIUS - RADIUS*c, WIDTH/2, -hd + RADIUS*c, lc)
    p[17] = geom.addPoint(LX + DX + RADIUS + RADIUS*c, WIDTH/2, -hd + RADIUS*c, lc)
    
    # -------------------------
    # INCIDENT
    # -------------------------
    p[18] = geom.addPoint(0, -WIDTH/2, -hd, lc)
    p[19] = geom.addPoint(0, WIDTH/2, -hd, lc)
    p[20] = geom.addPoint(0, -WIDTH/2, 0, lc)
    p[21] = geom.addPoint(0, WIDTH/2, 0, lc)
    
    # -------------------------
    # TRANSMISSION
    # -------------------------
    xT = 2*LX + 2*DX + 2*RADIUS
    
    p[22] = geom.addPoint(xT, -WIDTH/2, -hd, lc)
    p[23] = geom.addPoint(xT, WIDTH/2, -hd, lc)
    p[24] = geom.addPoint(xT, -WIDTH/2, 0, lc)
    p[25] = geom.addPoint(xT, WIDTH/2, 0, lc)

    ###############################
    # LINES
    # Store lines as a list
    l = [None]*42

    # -------------------------
    # ARCS
    # -------------------------
    l[0] = geom.addCircleArc(p[8],  p[9],  p[11])
    l[1] = geom.addCircleArc(p[11], p[9],  p[12])
    l[2] = geom.addCircleArc(p[12], p[9],  p[10])
    
    l[3] = geom.addCircleArc(p[13], p[14], p[16])
    l[4] = geom.addCircleArc(p[16], p[14], p[17])
    l[5] = geom.addCircleArc(p[17], p[14], p[15])
    
    # -------------------------
    # VERTICAL
    # -------------------------
    l[6]  = geom.addLine(p[20], p[18])
    l[7]  = geom.addLine(p[21], p[19])
    l[8]  = geom.addLine(p[4],  p[0])
    l[9]  = geom.addLine(p[6],  p[2])
    l[10] = geom.addLine(p[5],  p[1])
    l[11] = geom.addLine(p[7],  p[3])
    l[12] = geom.addLine(p[24], p[22])
    l[13] = geom.addLine(p[25], p[23])
    
    # -------------------------
    # HORIZONTAL Y
    # -------------------------
    l[14] = geom.addLine(p[18], p[19])
    l[15] = geom.addLine(p[20], p[21])
    l[16] = geom.addLine(p[8],  p[13])
    l[17] = geom.addLine(p[4],  p[6])
    l[18] = geom.addLine(p[10], p[15])
    l[19] = geom.addLine(p[5],  p[7])
    l[20] = geom.addLine(p[22], p[23])
    l[21] = geom.addLine(p[24], p[25])
    
    # -------------------------
    # HORIZONTAL X - INC & TRANS
    # -------------------------
    l[22] = geom.addLine(p[18], p[0])
    l[23] = geom.addLine(p[19], p[2])
    l[24] = geom.addLine(p[20], p[4])
    l[25] = geom.addLine(p[21], p[6])
    l[26] = geom.addLine(p[1],  p[22])
    l[27] = geom.addLine(p[3],  p[23])
    l[28] = geom.addLine(p[5],  p[24])
    l[29] = geom.addLine(p[7],  p[25])
    
    # -------------------------
    # TOP LINES IN ARC
    # -------------------------
    l[30] = geom.addLine(p[4], p[5])
    l[31] = geom.addLine(p[6], p[7])
    
    # -------------------------
    # DIAGONALS
    # -------------------------
    l[32] = geom.addLine(p[0], p[8])
    l[33] = geom.addLine(p[2], p[13])
    l[34] = geom.addLine(p[4], p[11])
    l[35] = geom.addLine(p[6], p[16])
    l[36] = geom.addLine(p[5], p[12])
    l[37] = geom.addLine(p[7], p[17])
    l[38] = geom.addLine(p[1], p[10])
    l[39] = geom.addLine(p[3], p[15])
    
    # -------------------------
    # HORIZONTAL Y - ARCS
    # -------------------------
    l[40] = geom.addLine(p[11], p[16])
    l[41] = geom.addLine(p[12], p[17])

    ###############################
    # CURVE LOOP
    # Store curve loops as a list of lists
    cl = [None]*18
    s  = [None]*18
    
    # ------------------------------------------------
    # ENTRANCE & OUTLET
    # ------------------------------------------------
    cl[0] = geom.addCurveLoop([ l[7], -l[14], -l[6],  l[15] ])
    cl[1] = geom.addCurveLoop([ l[12], l[20], -l[13], -l[21] ])
    
    s[0] = geom.addPlaneSurface([cl[0]])
    s[1] = geom.addPlaneSurface([cl[1]])
    
    # ------------------------------------------------
    # FREE SURFACE
    # ------------------------------------------------
    cl[2] = geom.addCurveLoop([ l[17], -l[25], -l[15],  l[24] ])
    cl[3] = geom.addCurveLoop([ l[19], -l[31], -l[17],  l[30] ])
    cl[4] = geom.addCurveLoop([ l[21], -l[29], -l[19],  l[28] ])
    
    s[2] = geom.addPlaneSurface([cl[2]])
    s[3] = geom.addPlaneSurface([cl[3]])
    s[4] = geom.addPlaneSurface([cl[4]])
    
    # ------------------------------------------------
    # RIGHT
    # ------------------------------------------------
    cl[5] = geom.addCurveLoop([ l[6],  l[22], -l[8],  -l[24] ])
    cl[6] = geom.addCurveLoop([ l[8],  l[32],  l[0],  -l[34] ])
    cl[7] = geom.addCurveLoop([ l[34], l[1],  -l[36], -l[30] ])
    cl[8] = geom.addCurveLoop([ l[36], l[2],  -l[38], -l[10] ])
    cl[9] = geom.addCurveLoop([ l[10], l[26], -l[12], -l[28] ])
    
    s[5] = geom.addPlaneSurface([cl[5]])
    s[6] = geom.addPlaneSurface([cl[6]])
    s[7] = geom.addPlaneSurface([cl[7]])
    s[8] = geom.addPlaneSurface([cl[8]])
    s[9] = geom.addPlaneSurface([cl[9]])

    # ------------------------------------------------
    # LEFT
    # ------------------------------------------------
    cl[10] = geom.addCurveLoop([ l[9],  -l[23], -l[7],  l[25] ])
    cl[11] = geom.addCurveLoop([ l[35], -l[3],  -l[33], -l[9] ])
    cl[12] = geom.addCurveLoop([ l[37], -l[4],  -l[35],  l[31] ])
    cl[13] = geom.addCurveLoop([ l[11],  l[39], -l[5],  -l[37] ])
    cl[14] = geom.addCurveLoop([ l[13], -l[27], -l[11],  l[29] ])
    
    s[10] = geom.addPlaneSurface([cl[10]])
    s[11] = geom.addPlaneSurface([cl[11]])
    s[12] = geom.addPlaneSurface([cl[12]])
    s[13] = geom.addPlaneSurface([cl[13]])
    s[14] = geom.addPlaneSurface([cl[14]])
    
    # ------------------------------------------------
    # ARCS
    # ------------------------------------------------
    cl[15] = geom.addCurveLoop([ l[16], l[3],  -l[40], -l[0] ])
    cl[16] = geom.addCurveLoop([ l[40], l[4],  -l[41], -l[1] ])
    cl[17] = geom.addCurveLoop([ l[41], l[5],  -l[18], -l[2] ])
    
    s[15] = geom.addSurfaceFilling([cl[15]])
    s[16] = geom.addSurfaceFilling([cl[16]])
    s[17] = geom.addSurfaceFilling([cl[17]])
    
    ###############################
    # TRANSFINITE
    # Transfinite surfaces
    for i in range(18):
        geom.mesh.setTransfiniteSurface(s[i])
    
    # Recombine surfaces (quad mesh)
    for i in range(18):
        geom.mesh.setRecombine(2, s[i])

    ###############################
    # TRANSFINITE
    # Transfinite lines
    # -------------------------
    # ARCS
    # -------------------------
    for i in [1, 4]:  # {2,5}
        geom.mesh.setTransfiniteCurve(l[i], ELEM, "Progression", 1)

    # TOP X LINES
    for i in [30, 31]:
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], ELEM, "Progression", 1)        
    
    # -------------------------
    # SHORT ARCS
    # -------------------------
    for i in [0, 2, 3, 5]:  # {1,3,4,6}
        geom.mesh.setTransfiniteCurve(l[i], ELEM_ARCS, "Progression", 1)
    
    # -------------------------
    # SQUARE BASE NEXT TO THE CYLINDER
    # -------------------------
    for i in range(32, 40):  # {33:40}
        geom.mesh.setTransfiniteCurve(l[i], 5, "Progression", 1)
    
    # -------------------------
    # REST OF THE DOMAIN X
    # -------------------------
    for i in range(22, 30):  # {23:30}
        geom.mesh.setTransfiniteCurve(l[i], ELEM * 3, "Progression", 1)
    
    # -------------------------
    # REST OF THE DOMAIN Y
    # -------------------------
    for i in range(14, 22):  # {15:22}
        geom.mesh.setTransfiniteCurve(l[i], 3, "Progression", 1)
    
    for i in range(40, 42):  # {41:42}
        geom.mesh.setTransfiniteCurve(l[i], 3, "Progression", 1)
    
    # -------------------------
    # HEIGHT OF THE DOMAIN
    # -------------------------
    for i in range(6, 14):  # {7:14}
        geom.mesh.setTransfiniteCurve(l[i], ELEM_ARCS, "Progression", 1)

    ###############################
    # PHYSICAL SURFACES 
    # Here we define a physical surface
    # ENTRANCE
    gmsh.model.addPhysicalGroup(2, [s[0]], 1)
    gmsh.model.setPhysicalName(2, 1, "ENTRANCE")
    
    # OUTLET
    gmsh.model.addPhysicalGroup(2, [s[1]], 2)
    gmsh.model.setPhysicalName(2, 2, "OUTLET")
    
    # POROUS_BOUNDARY {}
    gmsh.model.addPhysicalGroup(2, [], 3)
    gmsh.model.setPhysicalName(2, 3, "POROUS_BOUNDARY")
    
    # FREE_SURFACE {3:5}
    gmsh.model.addPhysicalGroup(2, s[2:5], 4)
    gmsh.model.setPhysicalName(2, 4, "FREE_SURFACE")
    
    # BOTTOM {}
    gmsh.model.addPhysicalGroup(2, [], 5)
    gmsh.model.setPhysicalName(2, 5, "BOTTOM")
    
    # IMPERMEABLE {16:18}
    gmsh.model.addPhysicalGroup(2, s[15:18], 6)
    gmsh.model.setPhysicalName(2, 6, "IMPERMEABLE")
    
    # FRONT {6:10}
    gmsh.model.addPhysicalGroup(2, s[5:10], 7)
    gmsh.model.setPhysicalName(2, 7, "FRONT")
    
    # BACK {11:15}
    gmsh.model.addPhysicalGroup(2, s[10:15], 8)
    gmsh.model.setPhysicalName(2, 8, "BACK")
    
    ###############################   

    # %% CREATE THE MESH:
    geom.synchronize()     

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'SEMI_CIRCULAR.msh'
    #filename = 'SEMI_CIRCULAR.m'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
