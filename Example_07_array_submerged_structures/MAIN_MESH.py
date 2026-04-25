import sys
import gmsh
import numpy as np
import os
import natsort
import math
from pathlib import Path

def MAIN_MESH(hd,FLUME_WIDTH,RADIO,HEIGHT,BATHY_X,LOC_STRUCT_X,LOC_STRUCT_Y):
    
    from OUTER_REGION import OUTER_REGION
       
    # LENGTH OF THE ELEMENTS IN GMSH
    lc  = 1.00
    
    # NUMBERS OF ELEMENTS FOR MESH
    ELEM_TOP = 6
      
    # BATHYMETRY 
    BATHY_Z = [-hd] * len(BATHY_X)
       
    path = os.getcwd()    
    
    #-----------------------------------#
    # SUB ROUTINES
    #-----------------------------------#
    OUTER_REGION(lc,FLUME_WIDTH,hd,RADIO,HEIGHT,LOC_STRUCT_X,LOC_STRUCT_Y,path,BATHY_X,BATHY_Z,ELEM_TOP)