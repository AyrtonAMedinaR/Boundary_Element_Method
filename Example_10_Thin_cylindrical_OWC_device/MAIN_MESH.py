import sys
import gmsh
import numpy as np
import os
import natsort
import math
from pathlib import Path

def MAIN_MESH(hd,FLUME_WIDTH,RADIO,DRAFT,LX,DX,LOC_STRUCT_X,LOC_STRUCT_Y):
    
    from REGION_INC import REGION_INC
    from REGION_OWC import REGION_OWC
    from REGION_OWC_CHAMBER import REGION_OWC_CHAMBER
    from REGION_TRA import REGION_TRA

    # LENGTH OF THE ELEMENTS IN GMSH
    lc  = 1.5
    
    # NUMBERS OF ELEMENTS FOR MESH
    ELEM_TOP = 6  # NUMBER OF ELEMENTS IN THE REEFS
    ELEM_X   = ELEM_TOP + 2

    path = os.getcwd()    
    
    #-----------------------------------#
    # SUB ROUTINES
    #-----------------------------------#
    REGION_INC(lc,ELEM_TOP,ELEM_X,hd,FLUME_WIDTH,RADIO,LX,DX,path)
    REGION_OWC(lc,ELEM_TOP,ELEM_X,hd,FLUME_WIDTH,RADIO,DRAFT,DX,LOC_STRUCT_X,LOC_STRUCT_Y,path)
    REGION_OWC_CHAMBER(lc,ELEM_TOP,ELEM_X,hd,FLUME_WIDTH,RADIO,DRAFT,DX,LOC_STRUCT_X,LOC_STRUCT_Y,path)
    REGION_TRA(lc,ELEM_TOP,ELEM_X,hd,FLUME_WIDTH,RADIO,LX,DX,path)

#os.chdir("../../")

