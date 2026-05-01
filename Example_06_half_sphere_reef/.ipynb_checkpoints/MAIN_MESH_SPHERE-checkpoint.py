import sys
import gmsh
import numpy as np
import os
import natsort
import math
from pathlib import Path

def MAIN_MESH_SPHERE(hd,FLUME_WIDTH,RADIO,BATHY_X,LOC_STRUCT_X,LOC_STRUCT_Y):
    
    from OUTER_REGION_SPHERE import OUTER_REGION_SPHERE
    from INNER_REGION_SPHERE import INNER_REGION_SPHERE
       
    # LENGTH OF THE ELEMENTS IN GMSH
    lc  = 0.10
    
    # NUMBERS OF ELEMENTS FOR MESH
    ELEM_TOP = 6
      
    # BATHYMETRY 
    BATHY_Z = [-hd] * len(BATHY_X)
       
    path = os.getcwd()    
    
    #-----------------------------------#
    # SUB ROUTINES
    #-----------------------------------#
    OUTER_REGION_SPHERE(lc,FLUME_WIDTH,hd,RADIO,LOC_STRUCT_X,LOC_STRUCT_Y,path,BATHY_X,BATHY_Z,ELEM_TOP)
    
    INNER_REGION_SPHERE(lc,hd,RADIO,LOC_STRUCT_X,LOC_STRUCT_Y,path,ELEM_TOP)