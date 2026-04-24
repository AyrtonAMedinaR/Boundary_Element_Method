import numpy as np

def BC_RADIATION(BC_RAD, ELEM_RAD, NE, Region):

    # Wave_height: Wave height
    # k: Wavenumber 
    # omega: Frequency
    # gravity: Gravitational accelaration
    # h: Water depth
    # POS: x, y and z location of each node
    # NCONEC: Number of nodes in a quad element 
    # ELEM_RAD: Elements belonging to the region where incident velocity potential is applied
    # KCONEC: Element connectivity array    
    # N: Total number of nodes in each region
    # NE: Total number of elements in each region
    # Region: Region where the BC is applied
    
    NE_prev = sum(NE[0:Region]) 

    for J in ELEM_RAD + NE_prev:

        BC_RAD[:, J] = 1.0 + 0*1j

    return BC_RAD
