import meshio
import numpy as np

def PROCESS_MESH_LINEAR(filename):
    # Reads a .msh file and extracts point coordinates and quad9 elements 
    # combined with their physical group data.
    
    # 1. Load the mesh file from the specified path
    mesh = meshio.read(filename)
    
    # 2. Extract point positions (coordinates)
    points = mesh.points
    
    # 3. Get the 'quad9' cell connectivity data
    quads = mesh.cells_dict["quad"]
    
    # 4. Get the physical group data associated with 'quad9' cells
    physical_data = mesh.cell_data_dict["gmsh:physical"]["quad"]
    
    # 5. Stack the physical data as a new column at the end of the connectivity matrix
    # Reshape(-1, 1) ensures the physical data is a column vector for hstack
    combined_quads = np.hstack([quads, physical_data.reshape(-1, 1)])
    
    return points, combined_quads
