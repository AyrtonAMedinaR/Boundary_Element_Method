def setRecombine(SV_params):
    import gmsh
    for localParam in SV_params:
        gmsh.model.mesh.setRecombine(2,localParam)

def INNER_REGION(lc,hd,RADIO,GAP,TRANS_X,TRANS_Y,path,ELEM_TOP):

    import gmsh
    from math import cos, sin, pi
    import numpy as np
    import os

    nStructures = int(len(TRANS_X))
    os.chdir(path)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    gmsh.model.add("modelo_2")

    # REDUCE CODE
    Sx = sin(45*pi/180)
    Cx = cos(45*pi/180)  

    ###############################
    # La siguiente línea permite abreviar el código:
    geom = gmsh.model.geo
    ###############################

    ###############################
    # --- BACK REEF
    ###############################
    # POINTS
    p = []

    for i in range(min(nStructures, len(TRANS_X))):
        p.append([
            # BASE
            geom.addPoint(TRANS_X[i],         TRANS_Y[i],                               GAP-hd,lc),
            geom.addPoint(( - RADIO)*Cx + TRANS_X[i], ( - RADIO)*Sx + TRANS_Y[i],       GAP-hd,lc),
            geom.addPoint(( - RADIO)*Sx + TRANS_X[i], (   RADIO)*Cx + TRANS_Y[i],       GAP-hd,lc),
            geom.addPoint((   RADIO)*Cx + TRANS_X[i], (   RADIO)*Sx + TRANS_Y[i],       GAP-hd,lc),
            geom.addPoint((   RADIO)*Sx + TRANS_X[i], ( - RADIO)*Cx + TRANS_Y[i],       GAP-hd,lc),
            #  TOP
            geom.addPoint(TRANS_X[i],         TRANS_Y[i],                               -hd,lc),
            geom.addPoint(( - RADIO)*Cx + TRANS_X[i], ( - RADIO)*Sx + TRANS_Y[i],       -hd,lc),
            geom.addPoint(( - RADIO)*Sx + TRANS_X[i], (   RADIO)*Cx + TRANS_Y[i],       -hd,lc),
            geom.addPoint((   RADIO)*Cx + TRANS_X[i], (   RADIO)*Sx + TRANS_Y[i],       -hd,lc),
            geom.addPoint((   RADIO)*Sx + TRANS_X[i], ( - RADIO)*Cx + TRANS_Y[i],      -hd,lc)
        ])

    # LINES
    ###############################
    # Store lines as a list
    l = []

    for i in range(min(nStructures, len(TRANS_X))):
        l.append([
            #
            geom.addCircleArc(p[i][6], p[i][5], p[i][9]),
            geom.addCircleArc(p[i][9], p[i][5], p[i][8]),
            geom.addCircleArc(p[i][8], p[i][5], p[i][7]),
            geom.addCircleArc(p[i][7], p[i][5], p[i][6]),
            #  ARCS BOTTOM
            geom.addCircleArc(p[i][1], p[i][0], p[i][4]),
            geom.addCircleArc(p[i][4], p[i][0], p[i][3]),
            geom.addCircleArc(p[i][3], p[i][0], p[i][2]),
            geom.addCircleArc(p[i][2], p[i][0], p[i][1]),
            #  SQUARE
            geom.addLine(p[i][1], p[i][4]),
            geom.addLine(p[i][4], p[i][3]),
            geom.addLine(p[i][3], p[i][2]),
            geom.addLine(p[i][2], p[i][1]),
            #  HEIGHT INTERIOR
            geom.addLine(p[i][1], p[i][6]),
            geom.addLine(p[i][4], p[i][9]),
            geom.addLine(p[i][3], p[i][8]),
            geom.addLine(p[i][2], p[i][7])
        ])

    # CURVE LOOP
    ###############################
    # Store curve loops as a list of lists
    cl = []

    for i in range(min(nStructures, len(TRANS_X))):
        cl.append([
            geom.addCurveLoop([l[i][8], l[i][9], l[i][10], l[i][11]]),
            geom.addCurveLoop([l[i][4], -l[i][8]]),
            geom.addCurveLoop([l[i][5], -l[i][9]]),
            geom.addCurveLoop([l[i][6], -l[i][10]]),
            geom.addCurveLoop([l[i][7], -l[i][11]]),
            #
            geom.addCurveLoop([l[i][12], l[i][0], -l[i][13], -l[i][4]]),
            geom.addCurveLoop([l[i][13], l[i][1], -l[i][14], -l[i][5]]),
            geom.addCurveLoop([l[i][14], l[i][2], -l[i][15], -l[i][6]]),
            geom.addCurveLoop([l[i][15], l[i][3], -l[i][12], -l[i][7]])
        ])

    # PLANE SURFACE
    ###############################
    # Store plane surfaces as a list of lists
    s = []

    for i in range(min(nStructures, len(TRANS_X))):
        s.append([
            geom.addPlaneSurface([cl[i][0]]),  
            geom.addPlaneSurface([cl[i][1]]),  
            geom.addPlaneSurface([cl[i][2]]),  
            geom.addPlaneSurface([cl[i][3]]),  
            geom.addPlaneSurface([cl[i][4]]),
            #
            geom.addSurfaceFilling([cl[i][5]]),
            geom.addSurfaceFilling([cl[i][6]]),
            geom.addSurfaceFilling([cl[i][7]]),
            geom.addSurfaceFilling([cl[i][8]])
        ])

    # TRANSFINITE CURVE
    ###############################
    #ELEM_TOP = 5  # NUMBER OF_B ELEMENTS IN THE REEFS

    for i in range(min(nStructures, len(TRANS_X))):  
        for j in range(8):  # Assuming you always have 20 curves per structure
            gmsh.model.geo.mesh.setTransfiniteCurve(l[i][j], ELEM_TOP + 2)

    #------- LAS ULTIMAS 5 SON DE LA BASE        
    for i in range(min(nStructures, len(TRANS_X))):  
        for j in range(12, 16):  # Indices for l31 to l35
            gmsh.model.geo.mesh.setTransfiniteCurve(l[i][j], ELEM_TOP)	
            
    for i in range(min(nStructures, len(TRANS_X))):  
        for j in range(8, 12):  # Indices for l21 to l30
            gmsh.model.geo.mesh.setTransfiniteCurve(l[i][j], ELEM_TOP)

    # TRANSFINITE SURFACE
    ###############################
    #-------
    for i in range(min(nStructures, len(TRANS_X))):
        for j in [0, 5, 6, 7, 8]:  # Iterate over all surfaces in structure i
            geom.mesh.setTransfiniteSurface(s[i][j])
    ###############################
    # --- BACK REEF
    ###############################            

    # PHYSICAL SURFACES
    ###############################   
    # Here we define a physical surface
    ps1 = gmsh.model.addPhysicalGroup(2, [ ], 1)
    gmsh.model.setPhysicalName(2, ps1, 'ENTRANCE')

    ps2 = gmsh.model.addPhysicalGroup(2, [ ], 2)
    gmsh.model.setPhysicalName(2, ps2, 'OUTLET')

    sum_1 = 0

    # BACK REEFS
    ###############################
    for i in range(min(nStructures, len(TRANS_X))):  
        if isinstance(s[i], list):
            num_surfaces = max(1, len(s[i]))
            sum_1 = 31 + i
            ps = gmsh.model.addPhysicalGroup(2, s[i][:num_surfaces], 31 + i)  # Use one ID per structure
            gmsh.model.setPhysicalName(2, ps, f'POROUS_BOUNDARY_{i+1}')
    ###############################       

    ps4 = gmsh.model.addPhysicalGroup(2,  [ ], 4)
    gmsh.model.setPhysicalName(2,  ps4, 'FREE_SURFACE')
   
    ps6 = gmsh.model.addPhysicalGroup(2, [ ], 6)
    gmsh.model.setPhysicalName(2, ps6, 'WALLS')

    ps7 = gmsh.model.addPhysicalGroup(2, [      ], 7)
    gmsh.model.setPhysicalName(2, ps7, 'FRONT_BOUNDARY')

    ps8 = gmsh.model.addPhysicalGroup(2, [      ], 8)
    gmsh.model.setPhysicalName(2, ps8, 'BACK_BOUNDARY')  

    ###############################   

    # %% CREATE THE MESH:
    geom.synchronize()

    # RECOMBINE
    ###############################          
    # BACK REEFS
    for i in range(min(nStructures, len(TRANS_X))):
        for j in range(len(s[i])):  # Iterate over all surfaces in structure i
            gmsh.model.mesh.setRecombine(2, s[i][j])       

    # SEE 2D ELEMENTS
    #gmsh.option.setNumber('Mesh.SurfaceFaces', 1)

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'INNER_REGION.msh'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
