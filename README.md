# 3D Boundary Element Method for Linear Water Waves

Python implementation of a 3D frequency–domain Boundary Element Method (BEM) to simulate linear wave interaction with submerged 
(impermeable and perforated) structures.

The code is designed for research in:

- Coastal protection  
- Artificial reefs  
- Porous and impermeable structures  
- Multi-domain wave interaction problems  
- Wave energy and OWC devices

## Reef ball simulation

<p align="center">
  <img src="docs/fig_wave_field.png" width="600">
</p>

The code computes hydrodynamic quantities such as:
- Velocity potential and fluxes
- Reflection, transmission and energy dissipation coefficients
- Wave forces

This repository is intended for research and educational use in coastal and ocean engineering. 
It includes a complete solver, utilities, and step-by-step tutorials.

---

# Main Features

• 3D frequency-domain Boundary Element Method  
• Multi-domain capability  
• Porous and impermeable boundaries  
• Constant, linear and quadratic boundary elements  
• Radiation and scattering problems  
• Hydrodynamic force computation  
• Energy flux and transmission analysis  
• Oblique wave incidence  
• Gmsh mesh workflow  
• Fully reproducible tutorials

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

# Repository Structure (Core Solver Modules)

BEM utilities for quadratic elements

```
BEM_UTILS/
```

General utilities used by the solver:

- Mesh processing
- Boundary identification
- Internal field evaluation
- Hydrodynamic forces
- Energy flux computation
- Free surface visualization

BEM utilities for constant elements

```
BEM_UTILS_CONSTANT/
```

BEM utilities for linear elements

```
BEM_UTILS_LINEAR/
```





