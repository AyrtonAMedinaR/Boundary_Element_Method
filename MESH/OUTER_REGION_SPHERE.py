def setRecombine(SV_params):
    import gmsh
    for localParam in SV_params:
        gmsh.model.mesh.setRecombine(2,localParam)

def OUTER_REGION_SPHERE(lc,FLUME_WIDTH,hd,RADIO,TRANS_X,TRANS_Y,path,BATHY_X,BATHY_Z,ELEM_TOP):

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
    # --- REEF POINTS
    ###############################    
    p1_2 = geom.addPoint(TRANS_X[0] + 0 ,           TRANS_Y[0] +  0  ,         -hd ,             lc)

    p2_2 = geom.addPoint(TRANS_X[0] + RADIO*Cx ,    TRANS_Y[0] + RADIO*Sx ,    -hd ,             lc)
    p3_2 = geom.addPoint(TRANS_X[0] - RADIO*Cx ,    TRANS_Y[0] + RADIO*Sx ,    -hd ,             lc)
    p4_2 = geom.addPoint(TRANS_X[0] + RADIO*Cx ,    TRANS_Y[0] - RADIO*Sx ,    -hd ,             lc)
    p5_2 = geom.addPoint(TRANS_X[0] - RADIO*Cx ,    TRANS_Y[0] - RADIO*Sx ,    -hd ,             lc)

    p6_2 = geom.addPoint(TRANS_X[0] - RADIO*Sx*Sx , TRANS_Y[0] + RADIO*Cx*Sx , -hd + RADIO*Sx ,  lc)
    p7_2 = geom.addPoint(TRANS_X[0] - RADIO*Sx*Sx , TRANS_Y[0] - RADIO*Cx*Sx , -hd + RADIO*Sx ,  lc)
    p8_2 = geom.addPoint(TRANS_X[0] + RADIO*Cx*Cx , TRANS_Y[0] + RADIO*Cx*Sx , -hd + RADIO*Sx ,  lc)
    p9_2 = geom.addPoint(TRANS_X[0] + RADIO*Cx*Cx , TRANS_Y[0] - RADIO*Cx*Sx , -hd + RADIO*Sx ,  lc)

    ###############################
    # 1 REEF LINES
    #  ARCS BOTTOM
    l1_2 = geom.addCircleArc(p4_2 , p1_2 , p5_2)
    l2_2 = geom.addCircleArc(p5_2 , p1_2 , p3_2)
    l3_2 = geom.addCircleArc(p3_2 , p1_2 , p2_2)
    l4_2 = geom.addCircleArc(p2_2 , p1_2 , p4_2)
    #  ARCS VERTICAL
    l5_2 = geom.addCircleArc(p4_2 , p1_2 , p9_2)
    l6_2 = geom.addCircleArc(p5_2 , p1_2 , p7_2)
    l7_2 = geom.addCircleArc(p3_2 , p1_2 , p6_2)
    l8_2 = geom.addCircleArc(p2_2 , p1_2 , p8_2)
    #  ARCS TOP
    l9_2 = geom.addCircleArc(p8_2 , p1_2 , p9_2)
    l10_2 = geom.addCircleArc(p6_2 , p1_2 , p8_2)
    l11_2 = geom.addCircleArc(p7_2 , p1_2 , p9_2)
    l12_2 = geom.addCircleArc(p6_2 , p1_2 , p7_2)

    ###############################
    # CURVE LOOP
    # CYLINDER 1
    cl1_2 = geom.addCurveLoop([l12_2, l11_2, -l9_2, -l10_2])
    cl2_2 = geom.addCurveLoop([l6_2, -l12_2, -l7_2, -l2_2])
    cl3_2 = geom.addCurveLoop([l7_2, l10_2, -l8_2, -l3_2])
    cl4_2 = geom.addCurveLoop([l8_2, l9_2, -l5_2, -l4_2])
    cl5_2 = geom.addCurveLoop([l5_2, -l11_2, -l6_2, -l1_2])

    s1_2  = geom.addSurfaceFilling([-cl1_2])
    s2_2  = geom.addSurfaceFilling([-cl2_2])
    s3_2  = geom.addSurfaceFilling([-cl3_2])
    s4_2  = geom.addSurfaceFilling([-cl4_2])    
    s5_2  = geom.addSurfaceFilling([-cl5_2]) 

    ##################
    # 1 REEF
    for tag in [l1_2, l2_2, l3_2, l4_2, l5_2, l6_2, l7_2, l8_2, l9_2, l10_2, l11_2, l12_2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP)

     ###############################          
    # 1 REEF
    geom.mesh.setTransfiniteSurface(s1_2)     
    geom.mesh.setTransfiniteSurface(s2_2)   
    geom.mesh.setTransfiniteSurface(s3_2)   
    geom.mesh.setTransfiniteSurface(s4_2)   
    geom.mesh.setTransfiniteSurface(s5_2)  
    
    ##################
    # BOX
    # VERTICAL LINES
    # for tag in [lp114, lp115, lp116, lp117]:
    #     gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP)

    # LATERAL LINES
    # for tag in [lp110, lp111, lp112, lp113]:
    #     gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP + 2)

    # HORIZONTAL LINES
    # for tag in [lp1, lp2, lp3, lp4]:
    #     gmsh.model.geo.mesh.setTransfiniteCurve(tag, 6*ELEM_TOP)

    ###############################        
    # %% CREATE THE MESH:
    geom.synchronize()

    ###############################            

    # PHYSICAL SURFACES
    ###############################   
    # Here we define a physical surface
    ps1 = gmsh.model.addPhysicalGroup(2, [ ps_e ], 1)
    gmsh.model.setPhysicalName(2, ps1, 'ENTRANCE')

    ps2 = gmsh.model.addPhysicalGroup(2, [ ps_o ], 2)
    gmsh.model.setPhysicalName(2, ps2, 'OUTLET')

    # CYL 1
    ps3 = gmsh.model.addPhysicalGroup(2,  [s1_2,s2_2,s3_2,s4_2,s5_2], 3)
    gmsh.model.setPhysicalName(2,  ps3, 'PERFORATED') 

    ps4 = gmsh.model.addPhysicalGroup(2,  [ ps_fs_1, ps_fs_2, ps_fs_3 ], 4)
    gmsh.model.setPhysicalName(2,  ps4, 'FREE_SURFACE')

    ps6 = gmsh.model.addPhysicalGroup(2, [ ], 6)
    gmsh.model.setPhysicalName(2, ps6, 'WALLS')

    ps7 = gmsh.model.addPhysicalGroup(2, [   ps_f_1, ps_f_2, ps_f_3   ], 7)
    gmsh.model.setPhysicalName(2, ps7, 'FRONT_BOUNDARY')

    ps8 = gmsh.model.addPhysicalGroup(2, [   ps_b_1, ps_b_2, ps_b_3   ], 8)
    gmsh.model.setPhysicalName(2, ps8, 'BACK_BOUNDARY')  
    ###############################        

    ###############################          
    # 1 REEF
    setRecombine([s1_2, s2_2, s3_2, s4_2, s5_2])    

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'OUTER_REGION_SPHERE.msh'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
