def setRecombine(SV_params):
    import gmsh
    for localParam in SV_params:
        gmsh.model.mesh.setRecombine(2,localParam)


def INNER_REGION_SPHERE(lc,hd,RADIO,TRANS_X,TRANS_Y,path,ELEM_TOP):

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
    # 1 REEF POINTS

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
    ###############################
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
    ###############################
    ###############################
    # CURVE LOOP
    # CYLINDER 1
    cl1_2 = geom.addCurveLoop([l12_2, l11_2, -l9_2, -l10_2])
    cl2_2 = geom.addCurveLoop([l6_2, -l12_2, -l7_2, -l2_2])
    cl3_2 = geom.addCurveLoop([l7_2, l10_2, -l8_2, -l3_2])
    cl4_2 = geom.addCurveLoop([l8_2, l9_2, -l5_2, -l4_2])
    cl5_2 = geom.addCurveLoop([l5_2, -l11_2, -l6_2, -l1_2])

    s1_2  = geom.addSurfaceFilling([cl1_2])
    s2_2  = geom.addSurfaceFilling([cl2_2])
    s3_2  = geom.addSurfaceFilling([cl3_2])
    s4_2  = geom.addSurfaceFilling([cl4_2])    
    s5_2  = geom.addSurfaceFilling([cl5_2]) 

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
    
    ###############################        
    # %% CREATE THE MESH:
    geom.synchronize()

    ############################### 
    # PHYSICAL NAME           
    # CYL 1
    ps3 = gmsh.model.addPhysicalGroup(2,  [s1_2,s2_2,s3_2,s4_2,s5_2], 3)
    gmsh.model.setPhysicalName(2,  ps3, 'PERFORATED') 

    ###############################          
    # 1 REEF
    setRecombine([s1_2, s2_2, s3_2, s4_2, s5_2])     

    ###############################          

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'INNER_REGION_SPHERE.msh'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
