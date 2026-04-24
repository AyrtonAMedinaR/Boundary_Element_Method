import sys
import gmsh
import numpy as np
import os
import natsort
import math
from pathlib import Path

def MAIN_MESH_SEMI_CIRC(hd,RADIUS,WIDTH,LX):
    
    from SEMI_CIRCULAR_OUTER import SEMI_CIRCULAR_OUTER
    from SEMI_CIRCULAR_INNER import SEMI_CIRCULAR_INNER

    # LENGTH OF THE ELEMENTS IN GMSH
    lc  = 0.50
    
    # NUMBERS OF ELEMENTS FOR MESH
    ELEM = 7;
    ELEM_ARCS = 7;
      
    DX = 0.10*hd;
     
    path = os.getcwd()    
    
    #-----------------------------------#
    # SUB ROUTINES
    #-----------------------------------#
    SEMI_CIRCULAR_OUTER(hd,RADIUS,WIDTH,LX,DX,lc,ELEM,ELEM_ARCS,path)
    SEMI_CIRCULAR_INNER(hd,RADIUS,WIDTH,LX,DX,lc,ELEM,ELEM_ARCS,path)