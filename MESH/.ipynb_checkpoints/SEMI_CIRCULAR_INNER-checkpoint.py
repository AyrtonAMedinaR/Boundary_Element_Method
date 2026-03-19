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

def SEMI_CIRCULAR_INNER(hd,RADIUS,WIDTH,LX,DX,lc,ELEM,ELEM_ARCS,path):

    os.chdir(path)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    gmsh.model.add("modelo_semi_circular_inner")

    ###############################
    # The following line allows you to abbreviate the code:
    geom = gmsh.model.geo
    ###############################

    ###############################
    # POINTS
    p = [None]*18
    c = math.cos(math.radians(45))

    # ARCS
    p[0]  = gmsh.model.geo.addPoint(LX + DX,            -WIDTH/2, -hd, lc)
    p[1]  = gmsh.model.geo.addPoint(LX + DX + RADIUS,   -WIDTH/2, -hd, lc)
    p[2]  = gmsh.model.geo.addPoint(LX + DX + 2*RADIUS, -WIDTH/2, -hd, lc)
    p[3]  = gmsh.model.geo.addPoint(LX + DX + RADIUS - RADIUS*c, -WIDTH/2, -hd + RADIUS*c, lc)
    p[4]  = gmsh.model.geo.addPoint(LX + DX + RADIUS + RADIUS*c, -WIDTH/2, -hd + RADIUS*c, lc)
    
    p[5]  = gmsh.model.geo.addPoint(LX + DX,            +WIDTH/2, -hd, lc)
    p[6]  = gmsh.model.geo.addPoint(LX + DX + RADIUS,   +WIDTH/2, -hd, lc)
    p[7]  = gmsh.model.geo.addPoint(LX + DX + 2*RADIUS, +WIDTH/2, -hd, lc)
    p[8]  = gmsh.model.geo.addPoint(LX + DX + RADIUS - RADIUS*c, +WIDTH/2, -hd + RADIUS*c, lc)
    p[9]  = gmsh.model.geo.addPoint(LX + DX + RADIUS + RADIUS*c, +WIDTH/2, -hd + RADIUS*c, lc)
    
    p[10] = gmsh.model.geo.addPoint(LX + DX + RADIUS/2,     -WIDTH/2, -hd, lc)
    p[11] = gmsh.model.geo.addPoint(LX + DX + 3*RADIUS/2,   -WIDTH/2, -hd, lc)
    p[12] = gmsh.model.geo.addPoint(LX + DX + RADIUS/2,     -WIDTH/2, -hd + RADIUS/2, lc)
    p[13] = gmsh.model.geo.addPoint(LX + DX + 3*RADIUS/2,   -WIDTH/2, -hd + RADIUS/2, lc)
    
    p[14] = gmsh.model.geo.addPoint(LX + DX + RADIUS/2,     +WIDTH/2, -hd, lc)
    p[15] = gmsh.model.geo.addPoint(LX + DX + 3*RADIUS/2,   +WIDTH/2, -hd, lc)
    p[16] = gmsh.model.geo.addPoint(LX + DX + RADIUS/2,     +WIDTH/2, -hd + RADIUS/2, lc)
    p[17] = gmsh.model.geo.addPoint(LX + DX + 3*RADIUS/2,   +WIDTH/2, -hd + RADIUS/2, lc)

    ###############################
    # LINES
    # Store lines as a list  
    l = [None]*26
    
    # ARCS
    l[0] = gmsh.model.geo.addCircleArc(p[0], p[1], p[3])
    l[1] = gmsh.model.geo.addCircleArc(p[3], p[1], p[4])
    l[2] = gmsh.model.geo.addCircleArc(p[4], p[1], p[2])
    l[3] = gmsh.model.geo.addCircleArc(p[5], p[6], p[8])
    l[4] = gmsh.model.geo.addCircleArc(p[8], p[6], p[9])
    l[5] = gmsh.model.geo.addCircleArc(p[9], p[6], p[7])
    
    # HORIZONTAL Y
    l[6] = gmsh.model.geo.addLine(p[0], p[5])
    l[7] = gmsh.model.geo.addLine(p[3], p[8])
    l[8] = gmsh.model.geo.addLine(p[4], p[9])
    l[9] = gmsh.model.geo.addLine(p[2], p[7])
    
    # FRONT DIAGONAL
    l[10] = gmsh.model.geo.addLine(p[0], p[10])
    l[11] = gmsh.model.geo.addLine(p[3], p[12])
    l[12] = gmsh.model.geo.addLine(p[4], p[13])
    l[13] = gmsh.model.geo.addLine(p[2], p[11])
    
    # FRONT HORIZONTAL
    l[14] = gmsh.model.geo.addLine(p[10], p[11])
    l[15] = gmsh.model.geo.addLine(p[12], p[13])
    
    # FRONT VERTICAL
    l[16] = gmsh.model.geo.addLine(p[12], p[10])
    l[17] = gmsh.model.geo.addLine(p[13], p[11])
    
    # BACK DIAGONAL
    l[18] = gmsh.model.geo.addLine(p[5], p[14])
    l[19] = gmsh.model.geo.addLine(p[8], p[16])
    l[20] = gmsh.model.geo.addLine(p[9], p[17])
    l[21] = gmsh.model.geo.addLine(p[7], p[15])
    
    # BACK HORIZONTAL
    l[22] = gmsh.model.geo.addLine(p[14], p[15])
    l[23] = gmsh.model.geo.addLine(p[16], p[17])
    
    # BACK VERTICAL
    l[24] = gmsh.model.geo.addLine(p[16], p[14])
    l[25] = gmsh.model.geo.addLine(p[17], p[15])

    ###############################
    # CURVE LOOP
    # Store curve loops as a list of lists
    cl = [None]*11
    s  = [None]*11
    
    # CURVED
    cl[0] = geom.addCurveLoop([ l[7], -l[3], -l[6],  l[0] ])
    cl[1] = geom.addCurveLoop([ l[8], -l[4], -l[7],  l[1] ])
    cl[2] = geom.addCurveLoop([ l[9], -l[5], -l[8],  l[2] ])
    
    s[0] = geom.addSurfaceFilling([cl[0]])
    s[1] = geom.addSurfaceFilling([cl[1]])
    s[2] = geom.addSurfaceFilling([cl[2]])
    
    # FRONT
    cl[3] = gmsh.model.geo.addCurveLoop([ l[10], -l[16], -l[11], -l[0] ])
    cl[4] = gmsh.model.geo.addCurveLoop([ l[15], -l[12], -l[1],  l[11] ])
    cl[5] = gmsh.model.geo.addCurveLoop([ l[12],  l[17], -l[13], -l[2] ])
    cl[6] = gmsh.model.geo.addCurveLoop([ l[16],  l[14], -l[17], -l[15] ])
    
    s[3] = gmsh.model.geo.addPlaneSurface([cl[3]])
    s[4] = gmsh.model.geo.addPlaneSurface([cl[4]])
    s[5] = gmsh.model.geo.addPlaneSurface([cl[5]])
    s[6] = gmsh.model.geo.addPlaneSurface([cl[6]])
    
    # BACK
    cl[7]  = gmsh.model.geo.addCurveLoop([ l[19],  l[24], -l[18],  l[3] ])
    cl[8]  = gmsh.model.geo.addCurveLoop([ l[20], -l[23], -l[19],  l[4] ])
    cl[9]  = gmsh.model.geo.addCurveLoop([ l[21], -l[25], -l[20],  l[5] ])
    cl[10] = gmsh.model.geo.addCurveLoop([ l[25], -l[22], -l[24],  l[23] ])
    
    s[7]  = gmsh.model.geo.addPlaneSurface([cl[7]])
    s[8]  = gmsh.model.geo.addPlaneSurface([cl[8]])
    s[9]  = gmsh.model.geo.addPlaneSurface([cl[9]])
    s[10] = gmsh.model.geo.addPlaneSurface([cl[10]])

    ###############################
    # TRANSFINITE
    # Transfinite surfaces
    for i in range(11):
        geom.mesh.setTransfiniteSurface(s[i])
    
    # Recombine surfaces (quad mesh)
    for i in range(11):
        geom.mesh.setRecombine(2, s[i])

    ###############################
    # TRANSFINITE
    # Transfinite lines
    # ARCS
    for i in [1, 4]:                 # {2,5}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], ELEM + 2, "Progression", 1)
    
    for i in range(14, 16):          # {15:16}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], ELEM + 2, "Progression",  1)
    
    for i in range(22, 24):          # {23:24}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], ELEM + 2, "Progression",  1)
    
    
    # SHORT ARCS
    for i in [0, 2, 3, 5]:           # {1,3,4,6}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], ELEM_ARCS, "Progression",  1)
    
    for i in range(16, 18):          # {17:18}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], ELEM_ARCS, "Progression",  1)
    
    for i in range(24, 26):          # {25:26}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], ELEM_ARCS, "Progression",  1)
    
    
    # SQUARE BASE NEXT TO THE CYLINDER
    for i in range(10, 14):          # {11:14}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], 5, "Progression",  1)
    
    for i in range(18, 22):          # {19:22}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], 5, "Progression",  1)
    
    
    # REST OF THE DOMAIN Y
    for i in range(6, 10):           # {7:10}
        gmsh.model.geo.mesh.setTransfiniteCurve(l[i], 3, "Progression",  1)



    ###############################
    # PHYSICAL SURFACES 
    # Here we define a physical surface

    gmsh.model.addPhysicalGroup(2, [], 1)
    gmsh.model.setPhysicalName(2, 1, "ENTRANCE")
    
    gmsh.model.addPhysicalGroup(2, [], 2)
    gmsh.model.setPhysicalName(2, 2, "OUTLET")
    
    gmsh.model.addPhysicalGroup(2, [s[0], s[1], s[2]], 3)
    gmsh.model.setPhysicalName(2, 3, "POROUS_BOUNDARY")
    
    gmsh.model.addPhysicalGroup(2, [], 4)
    gmsh.model.setPhysicalName(2, 4, "FREE_SURFACE")
    
    gmsh.model.addPhysicalGroup(2, [], 5)
    gmsh.model.setPhysicalName(2, 5, "BOTTOM")
    
    gmsh.model.addPhysicalGroup(2, [], 6)
    gmsh.model.setPhysicalName(2, 6, "IMPERMEABLE")
    
    gmsh.model.addPhysicalGroup(2, [s[3], s[4], s[5], s[6]], 7)
    gmsh.model.setPhysicalName(2, 7, "FRONT")
    
    gmsh.model.addPhysicalGroup(2, [s[7], s[8], s[9], s[10]], 8)
    gmsh.model.setPhysicalName(2, 8, "BACK")
    
    ###############################   

    # %% CREATE THE MESH:
    geom.synchronize()     

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'SEMI_CIRCULAR_INNER.msh'
    #filename = 'SEMI_CIRCULAR_INNER.m'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
