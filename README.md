# 3D Boundary Element Method for Linear Water Waves

Python implementation of a 3D frequency–domain Boundary Element Method (BEM) to simulate linear wave interaction with submerged  (impermeable and perforated) structures using quadratic (high-order) elements.

The code is designed for research in:
- Wave–structure interaction
- Coastal protection structures
- Artificial reefs
- OWC wave energy converter devices  

The code computes hydrodynamic quantities such as:
- Velocity potential and fluxes
- Reflection, transmission and energy dissipation coefficients
- Wave forces

This repository is intended for research and educational use in coastal and ocean engineering. 
It includes a complete solver, utilities, and step-by-step examples.

---

# Main Features

- 3D frequency-domain Boundary Element Method
- Higher-order elements discretisation for improved accuracy
- Constant and linear elements examples for comparison
- Multi-domain capability (for porous structures or large domains)
- Porous and impermeable boundary condition
- Radiation and scattering problems  
- Hydrodynamic force computation  
- Energy flux for reflection and transmission analysis  
- Oblique wave incidence  
- Mesh generation with Gmsh via .geo scripts or the Python API, with ready-to-use examples provided
- Fully reproducible examples

---

# Workflow

The workflow followed in all provided examples is:
- Define the fluid domain and generate the mesh using Gmsh
- Extract mesh and geometric characteristics
- Assemble and compute influence matrices
- Apply boundary conditions
- Post-process and analyse results

---

# Installation

## 1. Clone repository

```bash
git clone https://github.com/AyrtonAMedinaR/Boundary_Element_Method.git
cd Boundary_Element_Method
```
---

## 2. Create virtual environment

```
python -m venv bem_env
bem_env\Scripts\activate   # Windows
```

## 3. Install dependencies

```
pip install -r requirements.txt
```
---

# GMSH

Dependencies
Main packages:

numpy
scipy
matplotlib
gmsh (Python API)

You must install Gmsh separately:

```
https://gmsh.info
```

Make sure the Gmsh Python module works:

```
import gmsh
```

---

# Repository Structure (Core Solver Modules)

## 1. BEM utilities

General utilities used by the solver:

- Mesh processing
- Boundary identification
- Internal field evaluation
- Hydrodynamic forces
- Energy flux computation
- Free surface visualization

For quadratic elements

```
BEM_UTILS/
```

BEM utilities specific to constant elements

```
BEM_UTILS_CONSTANT/
```

BEM utilities specific to linear elements

```
BEM_UTILS_LINEAR/
```

## 2. Boundary Conditions

Implements all boundary conditions:

- Free surface
- Far-field radiation
- Scattering conditions
- Oblique wave incidence
- Interface coupling (multi-domain)
- Porous boundary formulation

Quadratic and linear element version:

```
BOUNDARY_CONDITIONS/
```

Constant-element version:

```
BOUNDARY_CONDITIONS_CONSTANT/
```

## 3. Influence Matrices (Core)

Main influence matrices:

- G and H matrices
- Numerical integration
- Fundamental solutions

The influence matrix functions implemented in this repository are based on the formulations presented in Domínguez (1993). 
The original implementation was provided in FORTRAN (for Elastodynamics applications), which has been both reformulated and simplified in Python for this work.

Domínguez, J. (1993). Boundary Elements in Dynamics: Computational Engineering. WIT Press.

```
MATRICES/
```

Constant elements

```
MATRICES_CONSTANT/
```

Linear elements

```
MATRICES_LINEAR/
```

---

# Examples

The repository includes complete step-by-step examples that serve as examples and validation cases.

## Example 01 — Empty Basin

Single and multi-domain wave propagation in an empty wave flume or basin.
Normal and oblique wave propagation can be considered.

```
Example_01_empty_wave_flume/
```

Outputs:

- Free surface elevation
- Reflection/transmission coefficients


## Example 02 — Submerged Semi-Circle (Impermeable)

Single-domain wave propagation over a submerged semi-circular obstacle.
Includes automated mesh generation using the Gmsh Python API.

```
Example_02_submerged_semi_circle_obstacle/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic forces 
- Free surface elevation
- Internal points calculation

## Example 03 — Multi-Domain Semi-Circle

Multi-domain wave propagation over a submerged semi-circular obstacle.

```
Example_03_multi_domain_submerged_semi_circle_obstacle/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic forces 
- Free surface elevation
- Internal points calculation


## Example 04 — Porous Semi-Circle (Multi-Domain)

Multi-domain wave propagation over a submerged porous semi-circular obstacle.

```
Example_04_multi_domain_porous_semi_circle_obstacle/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic forces 

## Example 05 — Porous Semi-Circle (Single Domain)

Single-domain wave propagation over a submerged porous semi-circular obstacle.
Includes automated mesh generation using the Gmsh Python API.

```
Example_05_porous_semi_circle_obstacle/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic forces 

## Example 06 — Half Sphere (porous and impermeable)

Single-domain wave propagation over a submerged porous/impermeable half-sphere obstacle.
Includes automated mesh generation using the Gmsh Python API.

```
Example_06_half_sphere_reef/
```

Outputs:

- Hydrodynamic forces 
- Free surface elevation

## Example 07 — Array of Impermeable Structures

Single-domain wave propagation over submerged impermeable cylindrical obstacles.
Includes automated mesh generation using the Gmsh Python API.

```
Example_07_array_submerged_structures/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic forces 
- Free surface elevation

## Example 08 — Array of Porous Structures

Multi-domain wave propagation over submerged porous cylindrical obstacles.
Includes automated mesh generation using the Gmsh Python API.

```
Example_08_array_porous_structures/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic forces 
- Free surface elevation

## Example 09 — Thick-Wall Cylindrical Oscillating Water Column (OWC)

Multi-domain wave interaction with a Thick-Wall Cylindrical Oscillating Water Column (OWC).
Includes automated mesh generation using the Gmsh Python API.

```
Example_09_Cylindrical_OWC_device/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic performance quantities


## Example 10 — Thin-Wall Cylindrical Oscillating Water Column (OWC)

Multi-domain wave interaction with a Thin-Wall Cylindrical Oscillating Water Column (OWC).
Includes automated mesh generation using the Gmsh Python API.

```
Example_10_Thin_cylindrical_OWC_device/
```

Outputs:

- Reflection/transmission coefficients
- Hydrodynamic performance quantities 
 
##  Constant and Linear Element Examples for empty wave flume water wave propagation

Comparison between discretization approaches.

```
Example_empty_wave_flume_CONSTANT_ELEMENTS/
```

```
Example_empty_wave_flume_LINEAR_ELEMENTS/
```

---

# How to Run a Simulation

1. Open any tutorial folder.
2. Launch the Jupyter notebook.
3. Run the notebook step-by-step.

Each tutorial contains:
- Mesh generation/loading
- Solver execution
- Post-processing and plots

---

# Citation

If you use this repository in academic work, please cite both the implementation and the original formulation source:

## This work
A. A. Medina Rodríguez (2026). *3D Boundary Element Method for Linear Water Waves*.

GitHub Repository: https://github.com/your-username/your-repo-name

## Theoretical foundation
Domínguez, J. (1993). *Boundary Elements in Dynamics: Computational Engineering* (Illustrated ed.). WIT Press.  
International Series on Computational Engineering / Topics in Engineering.  
ISBN: 978-1853122583.

---

# 🤝 Contributing

Any contributions you make are **greatly appreciated**.

---

# 📄 License

Distributed under the MIT License. See `LICENSE` for more information.

---

# 📬 Contact

**Ayrton Medina** [Email](mailto:ayrtonamedinar@gmail.com) | [GitHub](https://github.com/AyrtonAMedinaR) | [Google Scholar](https://scholar.google.com/citations?user=lUxvlWsAAAAJ) | [ResearchGate](https://www.researchgate.net/profile/Ayrton-Alfonso-Medina-Rodriguez?ev=hdr_xprf)

Project Link: [https://github.com/AyrtonAMedinaR/Boundary_Element_Method](https://github.com/AyrtonAMedinaR/Boundary_Element_Method)


