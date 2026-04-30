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

def REGION_OWC(lc,ELEM_TOP,ELEM_X,hd,FLUME_WIDTH,RADIO,DRAFT,DX,TRANS_X,TRANS_Y,path):

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
    p1_1 = geom.addPoint(TRANS_X[0],         TRANS_Y[0],       GAP-hd,lc)
    #  ARRIBA, CENTER
    p6_1 = geom.addPoint(TRANS_X[0],         TRANS_Y[0],                               0,lc)
    # MEDIO, EXTERIOR SHELL
    p11_1 = geom.addPoint(( - RADIO)*Cx + TRANS_X[0], ( - RADIO)*Sx + TRANS_Y[0],       GAP-hd,lc)
    p12_1 = geom.addPoint(( - RADIO)*Sx + TRANS_X[0], (   RADIO)*Cx + TRANS_Y[0],       GAP-hd,lc)
    p13_1 = geom.addPoint((   RADIO)*Cx + TRANS_X[0], (   RADIO)*Sx + TRANS_Y[0],       GAP-hd,lc)
    p14_1 = geom.addPoint((   RADIO)*Sx + TRANS_X[0], ( - RADIO)*Cx + TRANS_Y[0],       GAP-hd,lc)
    #  ARRIBA, EXTERIOR SHELL
    p15_1 = geom.addPoint(( - RADIO)*Cx + TRANS_X[0], ( - RADIO)*Sx + TRANS_Y[0],       0,lc)
    p16_1 = geom.addPoint(( - RADIO)*Sx + TRANS_X[0], (   RADIO)*Cx + TRANS_Y[0],       0,lc)
    p17_1 = geom.addPoint((   RADIO)*Cx + TRANS_X[0], (   RADIO)*Sx + TRANS_Y[0],       0,lc)
    p18_1 = geom.addPoint((   RADIO)*Sx + TRANS_X[0], ( - RADIO)*Cx + TRANS_Y[0],       0,lc)
    #  EXTERIOR SQUARE
    p19_1 = geom.addPoint(TRANS_X[0] - SQ_BOUND,  TRANS_Y[0] + SQ_BOUND,       0,lc)
    p20_1 = geom.addPoint(TRANS_X[0] - SQ_BOUND,  TRANS_Y[0] - SQ_BOUND,       0,lc)
    p21_1 = geom.addPoint(TRANS_X[0] + SQ_BOUND,  TRANS_Y[0] - SQ_BOUND,       0,lc)
    p22_1 = geom.addPoint(TRANS_X[0] + SQ_BOUND,  TRANS_Y[0] + SQ_BOUND,       0,lc)

    ###############################
    # 2 REEF POINTS
    p20_2 = geom.addPoint(TRANS_X[0] - SQ_BOUND,   FLUME_WIDTH/2,       0,lc)
    p21_2 = geom.addPoint(TRANS_X[0] + SQ_BOUND,   FLUME_WIDTH/2,       0,lc)

    ###############################
    # 4 REEF POINTS
    p19_4 = geom.addPoint(TRANS_X[0] - SQ_BOUND,  -FLUME_WIDTH/2,       0,lc)
    p22_4 = geom.addPoint(TRANS_X[0] + SQ_BOUND,  -FLUME_WIDTH/2,       0,lc)

    ###############################
    # BOX POINTS

    #  BOTTOM EXTERIOR SQUARE
    pb9_1  = geom.addPoint(TRANS_X[0] - SQ_BOUND,  TRANS_Y[0] + SQ_BOUND,        -hd,lc)
    pb10_1 = geom.addPoint(TRANS_X[0] - SQ_BOUND,  TRANS_Y[0] - SQ_BOUND,       -hd,lc)
    pb11_1 = geom.addPoint(TRANS_X[0] + SQ_BOUND,  TRANS_Y[0] - SQ_BOUND,       -hd,lc)
    pb12_1 = geom.addPoint(TRANS_X[0] + SQ_BOUND,  TRANS_Y[0] + SQ_BOUND,       -hd,lc)
	
    #  BOTTOM EXTERIOR SQUARE
    pb10_2 = geom.addPoint(TRANS_X[0] - SQ_BOUND,   FLUME_WIDTH/2,       -hd,lc)
    pb11_2 = geom.addPoint(TRANS_X[0] + SQ_BOUND,   FLUME_WIDTH/2,       -hd,lc)

    #  BOTTOM EXTERIOR SQUARE
    pb9_4  = geom.addPoint(TRANS_X[0] - SQ_BOUND,  -FLUME_WIDTH/2,        -hd,lc)
    pb12_4 = geom.addPoint(TRANS_X[0] + SQ_BOUND,  -FLUME_WIDTH/2,       -hd,lc)		 

    ###############################
    # LINES
    #  SQUARE
    l9_1 = geom.addLine(p11_1, p14_1)
    l10_1 = geom.addLine(p14_1, p13_1)
    l11_1 = geom.addLine(p13_1, p12_1)
    l12_1 = geom.addLine(p12_1, p11_1)

    #  ARCS EXTERIOR BOTTOM
    l21_1 = geom.addCircleArc(p11_1, p1_1, p12_1)
    l22_1 = geom.addCircleArc(p12_1, p1_1, p13_1)
    l23_1 = geom.addCircleArc(p13_1, p1_1, p14_1)
    l24_1 = geom.addCircleArc(p14_1, p1_1, p11_1)
    #  ARCS EXTERIOR TOP
    l25_1 = geom.addCircleArc(p15_1, p6_1, p16_1)
    l26_1 = geom.addCircleArc(p16_1, p6_1, p17_1)
    l27_1 = geom.addCircleArc(p17_1, p6_1, p18_1)
    l28_1 = geom.addCircleArc(p18_1, p6_1, p15_1)
    #  HEIGHT EXTERIOR
    l29_1 = geom.addLine(p15_1, p11_1)
    l30_1 = geom.addLine(p18_1, p14_1)
    l31_1 = geom.addLine(p17_1, p13_1)
    l32_1 = geom.addLine(p16_1, p12_1)
    #  SQUARE AROUND OWC
    l33_1 = geom.addLine(p19_1, p20_1)
    l34_1 = geom.addLine(p20_1, p21_1)
    l35_1 = geom.addLine(p21_1, p22_1)
    l36_1 = geom.addLine(p22_1, p19_1)
    l37_1 = geom.addLine(p16_1, p19_1)
    l38_1 = geom.addLine(p15_1, p20_1)
    l39_1 = geom.addLine(p18_1, p21_1)
    l40_1 = geom.addLine(p17_1, p22_1)

    ###############################
    # 2 REEF LINES

    l34_2 = geom.addLine(p20_2, p21_2)

    ###############################
    # 4 REEF LINES

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
    # 1 REEF CURVE LOOP
    cl1_1 = geom.addCurveLoop([l9_1, l10_1, l11_1, l12_1])
    cl2_1 = geom.addCurveLoop([l12_1, l21_1])
    cl3_1 = geom.addCurveLoop([l11_1, l22_1])
    cl4_1 = geom.addCurveLoop([l10_1, l23_1])
    cl5_1 = geom.addCurveLoop([l9_1,  l24_1])

    cl14_1 = geom.addCurveLoop([l29_1, -l24_1, -l30_1, l28_1])
    cl15_1 = geom.addCurveLoop([l30_1, -l23_1, -l31_1, l27_1])
    cl16_1 = geom.addCurveLoop([l31_1, -l22_1, -l32_1, l26_1])
    cl17_1 = geom.addCurveLoop([l32_1, -l21_1, -l29_1, l25_1])
    cl18_1 = geom.addCurveLoop([l33_1, -l38_1, l25_1, l37_1])
    cl19_1 = geom.addCurveLoop([l34_1, -l39_1, l28_1, l38_1])
    cl20_1 = geom.addCurveLoop([l35_1, -l40_1, l27_1, l39_1])
    cl21_1 = geom.addCurveLoop([l36_1, -l37_1, l26_1, l40_1])

    s1_1  = geom.addPlaneSurface([cl1_1]) 
    s2_1  = geom.addPlaneSurface([-cl2_1])
    s3_1  = geom.addPlaneSurface([-cl3_1])
    s4_1  = geom.addPlaneSurface([-cl4_1])
    s5_1  = geom.addPlaneSurface([-cl5_1])

    s14_1  = geom.addSurfaceFilling([-cl14_1])
    s15_1  = geom.addSurfaceFilling([-cl15_1])
    s16_1  = geom.addSurfaceFilling([-cl16_1])
    s17_1  = geom.addSurfaceFilling([-cl17_1]) 
    s18_1  = geom.addPlaneSurface([cl18_1])
    s19_1  = geom.addPlaneSurface([cl19_1])
    s20_1  = geom.addPlaneSurface([cl20_1])
    s21_1  = geom.addPlaneSurface([cl21_1]) 

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
    #PLANE SURFACE

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

    # ELEM_TOP = 5  # NUMBER OF ELEMENTS IN THE REEFS
    # ELEM_X   = 7

    # 1 REEF
    for tag in [l21_1, l22_1, l23_1, l24_1, l25_1, l26_1, l27_1, l28_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)

    # VERTICAL LINES    
    for tag in [l29_1, l30_1, l31_1, l32_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP)

    # SQUARE    
    for tag in [l9_1, l10_1, l11_1, l12_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_TOP)   
        
    # SQUARE AROUND OWC
    for tag in [l33_1, l34_1, l35_1, l36_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)  

    # ESPESOR    
    for tag in [l37_1, l38_1, l39_1, l40_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, 3)       

    for tag in [l34_2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)  

    for tag in [l36_4]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)  

    ##################
    # BOX

    # TOP OWC LINES FRONT        
    for tag in [lp33_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)
    # BOTTOM OWC LINES BACK        
    for tag in [lp35_1]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)

    # LATERALS BOX
    for tag in [ll_1,ll_2]:
        gmsh.model.geo.mesh.setTransfiniteCurve(tag, ELEM_X)        

    ###############################          
    # 1 REEF
    geom.mesh.setTransfiniteSurface(s1_1)     
 
    geom.mesh.setTransfiniteSurface(s14_1)   
    geom.mesh.setTransfiniteSurface(s15_1)  
    geom.mesh.setTransfiniteSurface(s16_1)  
    geom.mesh.setTransfiniteSurface(s17_1)  
    geom.mesh.setTransfiniteSurface(s18_1)  
    geom.mesh.setTransfiniteSurface(s19_1)  
    geom.mesh.setTransfiniteSurface(s20_1)  
    geom.mesh.setTransfiniteSurface(s21_1)           

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

    ps4 = gmsh.model.addPhysicalGroup(2,  [sb3_3,sb4_3,
                                           s18_1,s19_1,s20_1,s21_1], 4)
    gmsh.model.setPhysicalName(2,  ps4, 'FREE_SURFACE')

    ps5 = gmsh.model.addPhysicalGroup(2, [ ], 5)
    gmsh.model.setPhysicalName(2, ps5, 'BOTTOM')

    ps6 = gmsh.model.addPhysicalGroup(2, [s14_1,s15_1,s16_1,s17_1], 6)
    gmsh.model.setPhysicalName(2, ps6, 'IMPERMEABLE')    

    ps7 = gmsh.model.addPhysicalGroup(2, [sll_1], 7)
    gmsh.model.setPhysicalName(2, ps7, 'FRONT')

    ps8 = gmsh.model.addPhysicalGroup(2, [sll_2], 8)
    gmsh.model.setPhysicalName(2, ps8, 'BACK')   

    ps9 = gmsh.model.addPhysicalGroup(2,  [ ], 9)
    gmsh.model.setPhysicalName(2,  ps9, 'INTERNAL_FREE_SURFACE')    

    ps10 = gmsh.model.addPhysicalGroup(2,  [s1_1,s2_1,s3_1,s4_1,s5_1], 10)
    gmsh.model.setPhysicalName(2, ps10, 'MATCH_OWC')    

    ###############################          
    # 1 REEF
    setRecombine([s1_1, s2_1, s3_1, s4_1, s5_1,    s14_1, s15_1, s16_1, s17_1, s18_1, s19_1, s20_1, s21_1])

    setRecombine([sb5_1, sb6_1, sb7_1])   
    setRecombine([sb5_2, sb6_2, sb7_2])   
    setRecombine([sb3_3, sb4_3])   

    ## SEE 2D ELEMENTS
    ##gmsh.option.setNumber('Mesh.SurfaceFaces', 1)

    # SEE THE 1D NODES 
    gmsh.option.setNumber("Mesh.ElementOrder", 2)

    gmsh.model.mesh.generate(2)

    # SAVE THE MESH
    filename = 'REGION_OWC.msh'
    gmsh.write(filename)

    # OPEN THE INTERFACE TO SEE THE MESH IN GMSH
    gmsh.fltk.run()

    # END THE PROGRAM
    gmsh.finalize()
