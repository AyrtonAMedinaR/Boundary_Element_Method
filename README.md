# 3D Boundary Element Method for Linear Water Waves

Python implementation of a 3D frequency–domain Boundary Element Method (BEM) to simulate linear wave interaction with submerged 
(impermeable and perforated) structures using quadratic (high-order) elements.

The code is designed for research in:
- Wave–structure interaction
- Coastal protection structures
- Artificial reefs
- Wave energy and OWC devices
- Porous and impermeable structures  
- Multi-domain wave interaction problems  

The code computes hydrodynamic quantities such as:
- Velocity potential and fluxes
- Reflection, transmission and energy dissipation coefficients
- Wave forces

This repository is intended for research and educational use in coastal and ocean engineering. 
It includes a complete solver, utilities, and step-by-step tutorials.

---

# Main Features

- 3D frequency-domain Boundary Element Method  
- Multi-domain capability  
- Porous and impermeable boundaries  
- Constant, linear and quadratic boundary elements  
- Radiation and scattering problems  
- Hydrodynamic force computation  
- Energy flux for reflection and transmission analysis  
- Oblique wave incidence  
- Gmsh mesh workflow  
- Fully reproducible tutorials

---

# Theory Overview

The solver computes the velocity potential **φ** by solving the boundary integral equation for linear water waves.

It supports:

- Scattering problems  
- Radiation problems  
- Coupled multi-region domains  
- Porous media modelling via modified boundary conditions  

Outputs include:

- Reflection and transmission coefficients (CR, CT)  
- Free surface elevation  
- Hydrodynamic forces (Fx, Fy, Fz)  
- Energy flux across control surfaces  

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
-  Porous boundary formulation

Quadratic-element version:

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

Quadratic / higher order elements

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

# Tutorials

The repository includes complete step-by-step tutorials that serve as examples and validation cases.

## Tutorial 01 — Empty Basin

Single and multi-domain wave propagation in an empty wave flume or basin.
Normal and oblique wave propagation can be considered.

```
Tutorial_01_multi_domain_empty/
```

Outputs:

- Free surface validation
- Reflection/transmission coefficients


## Tutorial 02 — Submerged Semi-Circle (Impermeable)

```
Tutorial_02_submerged_semi_circle_impermeable/
```

Outputs:

- Free surface elevation
- Hydrodynamic forces

##  Tutorial 03 — Multi-Domain Semi-Circle

```
Tutorial_03_multi_domain_submerged_semi_circle_impermeable/
```

Demonstrates multi-region coupling.

##  Tutorial 04 — Porous Semi-Circle (Multi-Domain)

```
Tutorial_04_multi_domain_semi_circle_porous/
```

Demonstrates porous boundary modelling.

##  Tutorial 05 — Porous Semi-Circle (Single Domain)

```
Tutorial_05_semi_circle_porous/
```

Includes automated mesh generation.

##  Tutorial 06 — Half Sphere

```
Tutorial_06_half_sphere/
```

3D validation case:
- Impermeable sphere
- Porous sphere

##  Tutorial 07 — Array of Impermeable Structures

```
Tutorial_07_array_struct_impermeable/
```

Wave interaction with multiple bodies.

##  Tutorial 08 — Array of Porous Structures

```
Tutorial_08_array_struct_porous/
```

Demonstrates porous array behaviour.

##  Tutorial 09 — Cylindrical Oscillating Water Column (OWC)

```
Tutorial_09_Cylindrical_OWC/
```

Wave energy application:
- OWC performance metrics

##  Constant and Linear Element Examples

Comparison between discretization approaches.

```
Tutorial_CONSTANT_ELEMENTS/
```

```
Tutorial_LINEAR_ELEMENTS/
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

If you use this code in academic work, please cite:

Author: Ayrton A. Medina Rodríguez 3D Boundary Element Method for Linear Water Waves

---

# Contributing

Contributions and improvements are welcome.

Suggested contributions:
- Performance improvements
- Documentation updates
- New tutorial cases

---

# License

MIT

---

# Contact

For questions or collaboration:

Ayrton Medina (ayrtonamedinar@gmail.com)




