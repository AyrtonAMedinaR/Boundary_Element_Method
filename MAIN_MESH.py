import sys
import gmsh
import numpy as np
import os
import natsort
import math
from pathlib import Path

from OUTER_REGION import OUTER_REGION

from INNER_REGION import INNER_REGION
   
# LENGTH OF THE ELEMENTS IN GMSH
lc  = 0.75

# NUMBERS OF ELEMENTS FOR MESH
ELEM_TOP = 6

# WATER DEPTH
# BASED ON THE 3RD ROUND TEST AT OSU
# hd     = 3.0 # 0.697 0.767 0.830
# WAVE FLUME
FLUME_WIDTH = 10

# BACK REEF
# X Y
BR_LOC_X = []
BR_LOC_Y = []
with open("BR_X_Y.txt", "r") as file:
    for line in file:
        values = line.strip().split()  
        if len(values) == 2:  
            BR_LOC_X.append(float(values[0]))  
            BR_LOC_Y.append(float(values[1]))  

# BATHYMETRY 
BATHY_X = []
BATHY_Z = []
with open("BATHYMETRY.txt", "r") as file:
    for line in file:
        values = line.strip().split()  
        if len(values) == 2:  
            BATHY_X.append(float(values[0]))  
            BATHY_Z.append(float(values[1]))  


# CYLINDER DIMENSIONS hd,RADIO,GAP
hd     = abs(BATHY_Z[0])
RADIO  = 1*hd # (1/2)*hd
GAP    = (1 - 1/3)*hd


#-----------------------------------#
# FOLDERS ARE CREATED
#-----------------------------------#
# CASE = sys.argv[-1]
# # Directory
# CASE = CASE.zfill(2)
# directory = ("CASE_"+ CASE)
# #directory = "CASE_1"
    
# #Get CWD
# CurrentWD = os.getcwd()

# # Parent Directory path
# parent_dir = 'MESHES/'
# parent_path = os.path.join(CurrentWD, parent_dir)
# if os.path.exists(parent_path)==False: os.mkdir(parent_path)
# # Path
# path = os.path.join(CurrentWD, parent_dir, directory)
# # Create the directory
# # '/home / User / Documents'
# if os.path.exists(path)==False: 
#     os.mkdir(path)
#     #os.mkdir(path)
#     #print("Directory '% s' created" % directory)

path = os.getcwd()    

#-----------------------------------#
# SUB ROUTINES
#-----------------------------------#
OUTER_REGION(lc,FLUME_WIDTH,hd,RADIO,GAP,BR_LOC_X,BR_LOC_Y,path,BATHY_X,BATHY_Z,ELEM_TOP)

INNER_REGION(lc,hd,RADIO,GAP,BR_LOC_X,BR_LOC_Y,path,ELEM_TOP)

#os.chdir("../../")
