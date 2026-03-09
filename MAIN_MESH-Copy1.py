import sys
import gmsh
import numpy as np
import os
import natsort
import math
from pathlib import Path

def MAIN_MESH(hd,FLUME_WIDTH,RADIO,HEIGHT,BATHY_X,LOC_STRUCT_X,LOC_STRUCT_Y):
    
    from OUTER_REGION import OUTER_REGION
    from INNER_REGION import INNER_REGION
       
    # LENGTH OF THE ELEMENTS IN GMSH
    lc  = 0.75
    
    # NUMBERS OF ELEMENTS FOR MESH
    ELEM_TOP = 6
    
    # WAVE FLUME
    # FLUME_WIDTH = 10
    
    # BACK REEF
    # X Y  
    # BR_LOC_X = []
    # BR_LOC_Y = BR_LOC_X*0
    # with open("BR_X_Y.txt", "r") as file:
    #     for line in file:
    #         values = line.strip().split()  
    #         if len(values) == 2:  
    #             BR_LOC_X.append(float(values[0]))  
    #             BR_LOC_Y.append(float(values[1]))  
    
    # BATHYMETRY 
    BATHY_Z = [-hd] * len(BATHY_X)
    
    # BATHY_X = []
    # BATHY_Z = []
    # with open("BATHYMETRY.txt", "r") as file:
    #     for line in file:
    #         values = line.strip().split()  
    #         if len(values) == 2:  
    #             BATHY_X.append(float(values[0]))  
    #             BATHY_Z.append(float(values[1]))  
            
    # CYLINDER DIMENSIONS hd,RADIO,GAP
    # hd     = abs(BATHY_Z[0])
    # RADIO  = 1*hd # (1/2)*hd
    # GAP    = HEIGHT
    
    path = os.getcwd()    
    
    #-----------------------------------#
    # SUB ROUTINES
    #-----------------------------------#
    OUTER_REGION(lc,FLUME_WIDTH,hd,RADIO,HEIGHT,LOC_STRUCT_X,LOC_STRUCT_Y,path,BATHY_X,BATHY_Z,ELEM_TOP)
    
    INNER_REGION(lc,hd,RADIO,HEIGHT,LOC_STRUCT_X,LOC_STRUCT_Y,path,ELEM_TOP)