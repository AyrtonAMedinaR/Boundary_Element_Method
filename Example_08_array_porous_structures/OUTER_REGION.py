def setRecombine(SV_params):
    import gmsh
    for localParam in SV_params:
        gmsh.model.mesh.setRecombine(2,localParam)

def OUTER_REGION(lc,FLUME_WIDTH,hd,RADIO,GAP,TRANS_X,TRANS_Y,path,BATHY_X,BATHY_Z,ELEM_TOP):

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
    # --- BATHYMETRY
    ###############################
    point_front = []  # List to store point IDs
    point_back  = []  # List to store point IDs 

    # POINTS
    # BATHYMETRY 
    for i in range(len(BATHY_X)):
        tag_1 = i + 1 
        geom.addPoint(BATHY_X[i],   FLUME_WIDTH/2, BATHY_Z[i], lc )
        point_front.append(tag_1)            
    for j in range(len(BATHY_X)):
        tag_2 = len(BATHY_X) + j + 1 
        geom.addPoint(BATHY_X[j],  -FLUME_WIDTH/2, BATHY_Z[j], lc )
        point_back.append(tag_2) 
    # FREE SURFACE
    p_1 =  geom.addPoint(BATHY_X[0],    FLUME_WIDTH/2, 0, lc )
    p_2 =  geom.addPoint(BATHY_X[0],   -FLUME_WIDTH/2, 0, lc )     
    p_3 =  geom.addPoint(BATHY_X[-1],   FLUME_WIDTH/2, 0, lc )
    p_4 =  geom.addPoint(BATHY_X[-1],  -FLUME_WIDTH/2, 0, lc )  

    # FREE SURFACE POINTS TO DIVIDE INCIDENT, STURCTURES AND TRANSMITTED REGIONS
    p_5 =  geom.addPoint(BATHY_X[1],    FLUME_WIDTH/2, 0, lc )
    p_6 =  geom.addPoint(BATHY_X[1],   -FLUME_WIDTH/2, 0, lc )
    p_7 =  geom.addPoint(BATHY_X[-2],   FLUME_WIDTH/2, 0, lc )
    p_8 =  geom.addPoint(BATHY_X[-2],  -FLUME_WIDTH/2, 0, lc )           

    # LINES
    # BOTTOM 
    # f = front, b = bottom, e = entrance, o = outlet
    # d = down, u = up
    lext_f_d = geom.addPolyline(point_front[1:-1], 1) 
    lext_b_d = geom.addPolyline(point_back[1:-1],2)  
    lext_e_d = geom.addLine(point_front[0] , point_back[0]) 
    lext_o_d = geom.addLine(point_front[-1] , point_back[-1]) 
    # i = incident, t = transmitted
    lext_f_i_d = geom.addLine(point_front[0] , point_front[1]) 
    lext_b_i_d = geom.addLine(point_back[0] , point_back[1]) 
    lext_f_t_d = geom.addLine(point_front[-2] , point_front[-1]) 
    lext_b_t_d = geom.addLine(point_back[-2] , point_back[-1]) 
    # FREE SURFACE lines X direction
    lext_f_u_1 = geom.addLine(p_1 , p_5) 
    lext_f_u_2 = geom.addLine(p_5 , p_7) 
    lext_f_u_3 = geom.addLine(p_7 , p_3) 

    lext_b_u_1 = geom.addLine(p_2 , p_6)
    lext_b_u_2 = geom.addLine(p_6 , p_8)
    lext_b_u_3 = geom.addLine(p_8 , p_4)
    # FREE SURFACE lines Y direction
    lext_e_u = geom.addLine(p_1 , p_2)
    lext_o_u = geom.addLine(p_3 , p_4)
    # FREE SURFACE lines Y direction INTERMEDIATES
    lext_fs_u_1 = geom.addLine(p_5 , p_6)
    lext_fs_u_2 = geom.addLine(p_7 , p_8)
    # VERTICAL
    lext_1 = geom.addLine(point_front[0]  , p_1)
    lext_2 = geom.addLine(point_back[0]   , p_2) 
    lext_3 = geom.addLine(point_front[-1] , p_3)
    lext_4 = geom.addLine(point_back[-1]  , p_4)
    # VERTICAL INTERMEDIATE
    lext_5 = geom.addLine(point_front[1]  , p_5)
    lext_6 = geom.addLine(point_back[1]  , p_6)
    lext_7 = geom.addLine(point_front[-2]  , p_7)
    lext_8 = geom.addLine(point_back[-2]  , p_8)    

    # LOOPS FOR EXTERIOR BOUNDARIES
    cl_e = geom.addCurveLoop([ -lext_e_u,  -lext_1,  lext_e_d,  lext_2 ])
    cl_o = geom.addCurveLoop([  lext_o_u,  -lext_4, -lext_o_d,  lext_3 ])

    # FRONT WALL
    cl_f_1 = geom.addCurveLoop([  lext_f_u_1,  -lext_5, -lext_f_i_d,  lext_1 ])
    cl_f_2 = geom.addCurveLoop([  lext_f_u_2,  -lext_7, -lext_f_d,    lext_5 ])
    cl_f_3 = geom.addCurveLoop([  lext_f_u_3,  -lext_3, -lext_f_t_d,  lext_7 ])
    # BACK WALL
    cl_b_1 = geom.addCurveLoop([  lext_b_u_1,  -lext_6, -lext_b_i_d,  lext_2 ])
    cl_b_2 = geom.addCurveLoop([  lext_b_u_2,  -lext_8, -lext_b_d,    lext_6 ])
    cl_b_3 = geom.addCurveLoop([  lext_b_u_3,  -lext_4, -lext_b_t_d,  lext_8 ])
    # FREE SURFACE
    cl_fs_1 = geom.addCurveLoop([ -lext_f_u_1,  lext_e_u,     lext_b_u_1,  -lext_fs_u_1 ])
    cl_fs_2 = geom.addCurveLoop([ -lext_f_u_2,  lext_fs_u_1,  lext_b_u_2,  -lext_fs_u_2 ])
    cl_fs_3 = geom.addCurveLoop([ -lext_f_u_3,  lext_fs_u_2,  lext_b_u_3,  -lext_o_u ])

    # SURFACES
    ps_e  = geom.addPlaneSurface([ cl_e ])
    ps_o  = geom.addPlaneSurface([ cl_o ])
    # FRONT WALL
    ps_f_1  = geom.addPlaneSurface([ cl_f_1 ])
    ps_f_2  = geom.addPlaneSurface([ cl_f_2 ])
    ps_f_3  = geom.addPlaneSurface([ cl_f_3 ])
    # BACK WALL
    ps_b_1  = geom.addPlaneSurface([-cl_b_1 ])
    ps_b_2  = geom.addPlaneSurface([-cl_b_2 ])
    ps_b_3  = geom.addPlaneSurface([-cl_b_3 ])
    # FREE SURFACE
    ps_fs_1 = geom.addPlaneSurface([ cl_fs_1 ])
    ps_fs_2 = geom.addPlaneSurface([ cl_fs_2 ])
    ps_fs_3 = geom.addPlaneSurface([ cl_fs_3 ])

    # TRANSFINITE SURFACE
    geom.mesh.setTransfiniteSurface(ps_e)
    geom.mesh.setTransfiniteSurface(ps_o)
    # FRONT WALL
    geom.mesh.setTransfiniteSurface(ps_f_1)
    geom.mesh.setTransfiniteSurface(ps_f_2)
    geom.mesh.setTransfiniteSurface(ps_f_3)
    # BACK WALL
    geom.mesh.setTransfiniteSurface(ps_b_1)
    geom.mesh.setTransfiniteSurface(ps_b_2)
    geom.mesh.setTransfiniteSurface(ps_b_3)
    # FREE SURFACE
    geom.mesh.setTransfiniteSurface(ps_fs_1)
    geom.mesh.setTransfiniteSurface(ps_fs_2)
    geom.mesh.setTransfiniteSurface(ps_fs_3)
    ###############################
    # --- BATHYMETRY
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
            geom.addPlaneSurface([-cl[i][0]]),  
            geom.addPlaneSurface([-cl[i][1]]),  
            geom.addPlaneSurface([-cl[i][2]]),  
            geom.addPlaneSurface([-cl[i][3]]),  
            geom.addPlaneSurface([-cl[i][4]]),
            #
            geom.addSurfaceFilling([-cl[i][5]]),
            geom.addSurfaceFilling([-cl[i][6]]),
            geom.addSurfaceFilling([-cl[i][7]]),
            geom.addSurfaceFilling([-cl[i][8]])
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
    ps1 = gmsh.model.addPhysicalGroup(2, [ ps_e ], 1)
    gmsh.model.setPhysicalName(2, ps1, 'ENTRANCE')

    ps2 = gmsh.model.addPhysicalGroup(2, [ ps_o ], 2)
    gmsh.model.setPhysicalName(2, ps2, 'OUTLET')

    sum = 0
    for i in range(min(nStructures, len(TRANS_X))):
        sum = 31 + i
        ps = gmsh.model.addPhysicalGroup(2, s[i], 31 + i)  # Assign unique IDs dynamically
        gmsh.model.setPhysicalName(2, ps, f'POROUS_BOUNDARY_{i+1}')          

    ps4 = gmsh.model.addPhysicalGroup(2,  [ ps_fs_1, ps_fs_2, ps_fs_3 ], 4)
    gmsh.model.setPhysicalName(2,  ps4, 'FREE_SURFACE')

    ps6 = gmsh.model.addPhysicalGroup(2, [ ], 6)
    gmsh.model.setPhysicalName(2, ps6, 'WALLS')

    ps7 = gmsh.model.addPhysicalGroup(2, [   ps_f_1, ps_f_2, ps_f_3   ], 7)
    gmsh.model.setPhysicalName(2, ps7, 'FRONT_BOUNDARY')

    ps8 = gmsh.model.addPhysicalGroup(2, [   ps_b_1, ps_b_2, ps_b_3   ], 8)
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
    filename = 'OUTER_REGION.msh'
    #filename = 'OUTER_REGION.m'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()