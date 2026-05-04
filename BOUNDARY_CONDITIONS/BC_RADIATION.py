import numpy as np

def BC_RADIATION(BC_RAD, ELEM_RAD, NE, Region):

    # BC_RAD: Initizalize array
    # ELEM_RAD: Elements belonging to the region where the incident velocity potential is applied
    # NE: Total number of elements
    # Region: Region where the BC is applied
    
    NE_prev = sum(NE[0:Region]) 

    for J in ELEM_RAD + NE_prev:

        BC_RAD[:, J] = 1.0 + 0*1j

    return BC_RAD
