import sys
import gmsh
import numpy as np
import os
import natsort
from math import cos, sin, pi
from pathlib import Path

def setRecombine(SV_params):
    import gmsh
    for localParam in SV_params:
        gmsh.model.mesh.setRecombine(2,localParam)

def REGION_OWC_CHAMBER(lc,ELEM_TOP,ELEM_X,hd,FLUME_WIDTH,RADIO,DRAFT,DX,TRANS_X,TRANS_Y,path):

    os.chdir(path)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)

    gmsh.model.add("modelo_owc")

    ###############################
    # The following line allows you to abbreviate the code:
    geom = gmsh.model.geo
    ###############################

    Sx = sin(45*pi/180)
    Cx = cos(45*pi/180)

    POINT_INTER       = RADIO*Sx
    POINT_INTER_SHELL = RADIO*Sx
    SQ_BOUND          = RADIO + DX
    GAP               = hd - DRAFT         

    ###############################
    # 1 REEF POINTS

    #  MEDIO, CENTER
    p1_2 = geom.addPoint(TRANS_X[0],         TRANS_Y[0],                               0,lc)
    p2_2 = geom.addPoint(( - RADIO)*Cx + TRANS_X[0], ( - RADIO)*Sx + TRANS_Y[0],       0,lc)
    p3_2 = geom.addPoint(( - RADIO)*Sx + TRANS_X[0], (   RADIO)*Cx + TRANS_Y[0],       0,lc)
    p4_2 = geom.addPoint((   RADIO)*Cx + TRANS_X[0], (   RADIO)*Sx + TRANS_Y[0],       0,lc)
    p5_2 = geom.addPoint((   RADIO)*Sx + TRANS_X[0], ( - RADIO)*Cx + TRANS_Y[0],       0,lc)
    #  ARRIBA, CENTER
    p6_2 = geom.addPoint(TRANS_X[0],         TRANS_Y[0],                               GAP-hd,lc)
    p7_2 = geom.addPoint(( - RADIO)*Cx + TRANS_X[0], ( - RADIO)*Sx + TRANS_Y[0],       GAP-hd,lc)
    p8_2 = geom.addPoint(( - RADIO)*Sx + TRANS_X[0], (   RADIO)*Cx + TRANS_Y[0],       GAP-hd,lc)
    p9_2 = geom.addPoint((   RADIO)*Cx + TRANS_X[0], (   RADIO)*Sx + TRANS_Y[0],       GAP-hd,lc)
    p10_2 = geom.addPoint((   RADIO)*Sx + TRANS_X[0], ( - RADIO)*Cx + TRANS_Y[0],      GAP-hd,lc)

    ###############################
    ###############################
    ###############################
    # 1 REEF LINES
    l1_2 = geom.addCircleArc(p7_2, p6_2, p10_2)
    l2_2 = geom.addCircleArc(p10_2, p6_2, p9_2)
    l3_2 = geom.addCircleArc(p9_2, p6_2, p8_2)
    l4_2 = geom.addCircleArc(p8_2, p6_2, p7_2)
    #  ARCS BOTTOM
    l5_2 = geom.addCircleArc(p2_2, p1_2, p5_2)
    l6_2 = geom.addCircleArc(p5_2, p1_2, p4_2)
    l7_2 = geom.addCircleArc(p4_2, p1_2, p3_2)
    l8_2 = geom.addCircleArc(p3_2, p1_2, p2_2)
    #  SQUARE
    l9_2 = geom.addLine(p2_2,  p5_2)
    l10_2 = geom.addLine(p5_2, p4_2)
    l11_2 = geom.addLine(p4_2, p3_2)
    l12_2 = geom.addLine(p3_2, p2_2)
    #  HEIGHT INTERIOR
    l13_2 = geom.addLine(p2_2, p7_2)
    l14_2 = geom.addLine(p5_2, p10_2)
    l15_2 = geom.addLine(p4_2, p9_2)
    l16_2 = geom.addLine(p3_2, p8_2)     
    #  SQUARE GAP
    l17_2 = geom.addLine(p7_2, p8_2)
    l18_2 = geom.addLine(p8_2, p9_2)
    l19_2 = geom.addLine(p9_2, p10_2)
    l20_2 = geom.addLine(p10_2, p7_2)

    ###############################
    ###############################
    ###############################
    # CURVE LOOP
    # CYLINDER 1
    cl1_2 = geom.addCurveLoop([l9_2, l10_2, l11_2, l12_2])
    cl2_2 = geom.addCurveLoop([l5_2, -l9_2])
    cl3_2 = geom.addCurveLoop([l6_2, -l10_2])
    cl4_2 = geom.addCurveLoop([l7_2, -l11_2])
    cl5_2 = geom.addCurveLoop([l8_2, -l12_2])

    cl6_2 = geom.addCurveLoop([l13_2, l1_2, -l14_2, -l5_2])
    cl7_2 = geom.addCurveLoop([l14_2, l2_2, -l15_2, -l6_2])
    cl8_2 = geom.addCurveLoop([l15_2, l3_2, -l16_2, -l7_2])
    cl9_2 = geom.addCurveLoop([l16_2, l4_2, -l13_2, -l8_2])

    cl10_2 = geom.addCurveLoop([l17_2, l18_2, l19_2, l20_2])
    cl11_2 = geom.addCurveLoop([l1_2, l20_2])
    cl12_2 = geom.addCurveLoop([l2_2, l19_2])
    cl13_2 = geom.addCurveLoop([l3_2, l18_2])
    cl14_2 = geom.addCurveLoop([l4_2, l17_2])



    s1_2  = geom.addPlaneSurface([cl1_2]) 
    s2_2  = geom.addPlaneSurface([cl2_2])
    s3_2  = geom.addPlaneSurface([cl3_2])
    s4_2  = geom.addPlaneSurface([cl4_2])
    s5_2  = geom.addPlaneSurface([cl5_2])    

    s6_2  = geom.addSurfaceFilling([cl6_2])
    s7_2  = geom.addSurfaceFilling([cl7_2])
    s8_2  = geom.addSurfaceFilling([cl8_2])
    s9_2  = geom.addSurfaceFilling([cl9_2])

    s10_2  = geom.addPlaneSurface([cl10_2]) 
    s11_2  = geom.addPlaneSurface([-cl11_2]) 
    s12_2  = geom.addPlaneSurface([-cl12_2]) 
    s13_2  = geom.addPlaneSurface([-cl13_2]) 
    s14_2  = geom.addPlaneSurface([-cl14_2]) 

    ##################
    # 1 REEF
    for tag in [l1_2, l2_2, l3_2, l4_2, l5_2, l6_2, l7_2, l8_2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP + 2)
    for tag in [l13_2, l14_2, l15_2, l16_2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP)
    for tag in [l9_2, l10_2, l11_2, l12_2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP)    

    for tag in [l17_2, l18_2, l19_2, l20_2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP)

    ###############################          
    # 1 REEF
    geom.mesh.setTransfiniteSurface(s1_2)  

    geom.mesh.setTransfiniteSurface(s6_2)   
    geom.mesh.setTransfiniteSurface(s7_2)   
    geom.mesh.setTransfiniteSurface(s8_2)   
    geom.mesh.setTransfiniteSurface(s9_2)

    geom.mesh.setTransfiniteSurface(s10_2)   

    ###############################        
    # %% CREATE THE MESH:
    geom.synchronize()

    ############################### 
    # PHYSICAL NAME           

    ps6 = gmsh.model.addPhysicalGroup(2, [s6_2,s7_2,s8_2,s9_2], 6)
    gmsh.model.setPhysicalName(2, ps6, 'IMPERMEABLE')  

    # CYL 1
    ps9 = gmsh.model.addPhysicalGroup(2,  [s1_2,s2_2,s3_2,s4_2,s5_2], 9)
    gmsh.model.setPhysicalName(2,  ps9, 'INTERNAL_FREE_SURFACE')
    
    ps10 = gmsh.model.addPhysicalGroup(2,  [s10_2,s11_2,s12_2,s13_2,s14_2], 10)
    gmsh.model.setPhysicalName(2, ps10, 'MATCH_OWC')

    ###############################          
    # 1 REEF
    setRecombine([s1_2, s2_2, s3_2, s4_2, s5_2, s6_2, s7_2, s8_2, s9_2])
    setRecombine([s10_2, s11_2, s12_2, s13_2, s14_2])

    ###############################          

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'REGION_OWC_CHAMBER.msh'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
