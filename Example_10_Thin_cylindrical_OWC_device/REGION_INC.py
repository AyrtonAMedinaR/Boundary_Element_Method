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

def REGION_INC(lc,ELEM_TOP,ELEM_X,hd,FLUME_WIDTH,RADIO,LENGTH_X,DX,path):

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

    ###############################
    # POINTS, TOP

    #  EXTERIOR SQUARE
    p19_1 = geom.addPoint(0,         + SQ_BOUND,       0,lc)
    p20_1 = geom.addPoint(0,         - SQ_BOUND,       0,lc)
    p21_1 = geom.addPoint(LENGTH_X,  - SQ_BOUND,       0,lc)
    p22_1 = geom.addPoint(LENGTH_X,  + SQ_BOUND,       0,lc)

    p20_2 = geom.addPoint(0,         FLUME_WIDTH/2,       0,lc)
    p21_2 = geom.addPoint(LENGTH_X,  FLUME_WIDTH/2,       0,lc)

    p19_4 = geom.addPoint(0,        -FLUME_WIDTH/2,       0,lc)
    p22_4 = geom.addPoint(LENGTH_X, -FLUME_WIDTH/2,       0,lc)

    ###############################
    # POINTS, BOTTOM

    #  EXTERIOR SQUARE
    pb9_1  = geom.addPoint(0,         + SQ_BOUND,       -hd,lc)
    pb10_1 = geom.addPoint(0,         - SQ_BOUND,       -hd,lc)
    pb11_1 = geom.addPoint(LENGTH_X,  - SQ_BOUND,       -hd,lc)
    pb12_1 = geom.addPoint(LENGTH_X,  + SQ_BOUND,       -hd,lc)

    pb10_2 = geom.addPoint(0,         FLUME_WIDTH/2,       -hd,lc)
    pb11_2 = geom.addPoint(LENGTH_X,  FLUME_WIDTH/2,       -hd,lc)

    pb9_4  = geom.addPoint(0,        -FLUME_WIDTH/2,       -hd,lc)
    pb12_4 = geom.addPoint(LENGTH_X, -FLUME_WIDTH/2,       -hd,lc)   

    ###############################
    # LINES
    
    #  SQUARE AROUND OWC
    l33_1 = geom.addLine(p19_1, p20_1)
    l34_1 = geom.addLine(p20_1, p21_1)
    l35_1 = geom.addLine(p21_1, p22_1)
    l36_1 = geom.addLine(p22_1, p19_1)

    l34_2 = geom.addLine(p20_2, p21_2)

    l36_4 = geom.addLine(p22_4, p19_4)

    ###############################
    # BOX LINES

    # OWC bottom lines (the 2nd and 3rd were deleted because the bottom is not meshed)
    lp33_1 = geom.addLine(pb9_1, pb10_1)

    lp35_1 = geom.addLine(pb11_1, pb12_1)


    # OWC 2 AND 1
    lp_1_2 = geom.addLine(p20_2, p19_1)
    lp_2_2 = geom.addLine(p21_2, p22_1)    

    # OWC 1 AND 4
    lp_1_1 = geom.addLine(p20_1, p19_4)
    lp_2_1 = geom.addLine(p21_1, p22_4)  

    # OWC 2 AND 1
    lpp_1_2 = geom.addLine(pb10_2, pb9_1)
    lpp_2_2 = geom.addLine(pb11_2, pb12_1)            

    # OWC 1 AND 4
    lpp_1_1 = geom.addLine(pb10_1, pb9_4)
    lpp_2_1 = geom.addLine(pb11_1, pb12_4)  

    # VERTICAL LINES
    # OWC 1
    lv_1_1 = geom.addLine(p19_1, pb9_1)
    lv_2_1 = geom.addLine(p20_1, pb10_1)
    lv_3_1 = geom.addLine(p21_1, pb11_1)
    lv_4_1 = geom.addLine(p22_1, pb12_1)

    # OWC 2
    lv_2_2 = geom.addLine(p20_2, pb10_2)
    lv_3_2 = geom.addLine(p21_2, pb11_2)

    # OWC 4
    lv_1_4 = geom.addLine(p19_4, pb9_4)
    lv_4_4 = geom.addLine(p22_4, pb12_4)

    # LATERALES INFERIORES
    ll_1 = geom.addLine(pb9_4, pb12_4)
    ll_2 = geom.addLine(pb10_2, pb11_2)   

    ###############################
    # CURVE LOOP
    cl1_1 = geom.addCurveLoop([l33_1, l34_1, l35_1, l36_1])

    s1_1  = geom.addPlaneSurface([cl1_1]) 

    ###############################
    # BOX CURVE LOOP

    # FRONT
    clb5_1  = geom.addCurveLoop([lv_2_2, lpp_1_2, -lv_1_1, -lp_1_2]) 
    clb6_1  = geom.addCurveLoop([lv_1_1, lp33_1,  -lv_2_1, -l33_1]) 
    clb7_1  = geom.addCurveLoop([lv_2_1, lpp_1_1, -lv_1_4, -lp_1_1]) 

    # BACK
    clb5_2  = geom.addCurveLoop([lv_3_2,  lpp_2_2, -lv_4_1, -lp_2_2]) 
    clb6_2  = geom.addCurveLoop([lv_4_1, -lp35_1,  -lv_3_1,  l35_1]) 
    clb7_2  = geom.addCurveLoop([lv_3_1,  lpp_2_1, -lv_4_4, -lp_2_1]) 

    # TOP
    clb3_3  = geom.addCurveLoop([-l34_2,  lp_1_2, -l36_1, -lp_2_2])
    clb4_3  = geom.addCurveLoop([-l34_1,  lp_1_1, -l36_4, -lp_2_1])

    # LATERALS
    cll_1 = geom.addCurveLoop([l36_4, lv_1_4, ll_1, -lv_4_4])
    cll_2 = geom.addCurveLoop([l34_2, lv_3_2,-ll_2, -lv_2_2])    

    ###############################
    # PLANE SURFACE

    sb5_1 = geom.addPlaneSurface([clb5_1])
    sb6_1 = geom.addPlaneSurface([clb6_1])
    sb7_1 = geom.addPlaneSurface([clb7_1])

    sb5_2 = geom.addPlaneSurface([-clb5_2])
    sb6_2 = geom.addPlaneSurface([-clb6_2])
    sb7_2 = geom.addPlaneSurface([-clb7_2]) 

    sb3_3 = geom.addPlaneSurface([clb3_3])
    sb4_3 = geom.addPlaneSurface([clb4_3])  

    sll_1 = geom.addPlaneSurface([cll_1])  
    sll_2 = geom.addPlaneSurface([cll_2])           

    ###############################
    # TRANSFINITE CURVE
    for tag in [l33_1, l35_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)      

    ##################
    # BOX

    # TOP OWC LINES FRONT Y       
    for tag in [lp33_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)
    # BOTTOM OWC LINES BACK Y       
    for tag in [lp35_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)

    ###############################          
    # TRANSFINITE SURFACE
    geom.mesh.setTransfiniteSurface(s1_1)           

    # FRONT
    geom.mesh.setTransfiniteSurface(sb5_1) 
    geom.mesh.setTransfiniteSurface(sb6_1) 
    geom.mesh.setTransfiniteSurface(sb7_1)     

    # BACK
    geom.mesh.setTransfiniteSurface(sb5_2) 
    geom.mesh.setTransfiniteSurface(sb6_2) 
    geom.mesh.setTransfiniteSurface(sb7_2)    

    # TOP
    geom.mesh.setTransfiniteSurface(sb3_3)
    geom.mesh.setTransfiniteSurface(sb4_3)  

    # LATERALS
    geom.mesh.setTransfiniteSurface(sll_1) 
    geom.mesh.setTransfiniteSurface(sll_2)          

    # %% CREATE THE MESH:
    geom.synchronize()

    ###############################   

    # Here we define a physical surface
    ps1 = gmsh.model.addPhysicalGroup(2, [sb5_1,sb6_1,sb7_1], 1)
    gmsh.model.setPhysicalName(2, ps1, 'ENTRANCE')

    ps2 = gmsh.model.addPhysicalGroup(2, [sb5_2,sb6_2,sb7_2], 2)
    gmsh.model.setPhysicalName(2, ps2, 'OUTLET')

    # OWC 1 - 5
    ps4 = gmsh.model.addPhysicalGroup(2,  [s1_1,   sb3_3,sb4_3], 4)
    gmsh.model.setPhysicalName(2,  ps4, 'FREE_SURFACE')

    ps5 = gmsh.model.addPhysicalGroup(2, [ ], 5)
    gmsh.model.setPhysicalName(2, ps5, 'BOTTOM')

    # OWC 1 - 5
    ps6 = gmsh.model.addPhysicalGroup(2, [ ], 6)
    gmsh.model.setPhysicalName(2, ps6, 'IMPERMEABLE')

    ps7 = gmsh.model.addPhysicalGroup(2, [sll_1], 7)
    gmsh.model.setPhysicalName(2, ps7, 'FRONT')

    ps8 = gmsh.model.addPhysicalGroup(2, [sll_2], 8)
    gmsh.model.setPhysicalName(2, ps8, 'BACK')  

    ###############################          
    # RECOMBINE
    setRecombine([s1_1])

    setRecombine([sb5_1, sb6_1, sb7_1])   
    setRecombine([sb5_2, sb6_2, sb7_2])   
    setRecombine([sb3_3, sb4_3])   

    ## SEE 2D ELEMENTS
    ##gmsh.option.setNumber('Mesh.SurfaceFaces', 1)

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'REGION_INC.msh'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
