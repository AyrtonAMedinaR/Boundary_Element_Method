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

<p align="center">
  <img width="600" height="552" alt="RANDOM_LAYOUT" src="https://github.com/user-attachments/assets/c79ad3ef-f7a6-41cc-93c8-5a75ff619db7" />
  <br>
  <em>Wave interaction with a random layout of submerged perforated cylinders simulated using the 3D BEM (Medina Rodriguez et al., 2026).</em>
</p>

## Publications

The following peer-reviewed publications use or are directly related to this 3D Boundary Element Method.

**Medina Rodríguez, A.A., Huang, Z., R3D Consortium (2026)**  
Wave Scattering and Energy Reduction by a Patch of Idealized Coral-Growing Units.  
In: *Coastal Dynamics 2025*, Coastal Research Library, Vol. 41, Springer.  
https://doi.org/10.1007/978-3-032-15473-6_85

**Medina Rodríguez, A.A., Trivedi, K., Koley, S., Oderiz Martinez, I., Mendoza, E., Posada Vanegas, G., Silva, R. (2023)**  
Improved hydrodynamic performance of an OWC device based on a Helmholtz resonator.  
*Energy*, Volume 273, 127299.  
https://doi.org/10.1016/j.energy.2023.127299

**Medina Rodríguez, A.A., Silva Casarín, R., Blanco Ilzarbe, J.M. (2022)**  
A 3D boundary element method for analysing the hydrodynamic performance of a land-fixed oscillating water column device.  
*Engineering Analysis with Boundary Elements*, Volume 138, pp. 407-422.  
https://doi.org/10.1016/j.enganabound.2022.02.014

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

# Theory

The Cartesian coordinate system is adopted and the origin is assumed at the undisturbed free surface water level. Waves propagate and approach the structures from the left in the $x$ direction, in water of depth $DEPTH$. The water flow is considered incompressible and irrotational, and the effect of surface tension is neglected. Linear water wave theory and linear hydrodynamic interactions are investigated under the assumption that the wave amplitude is much smaller than the wavelength.  A simple harmonic flow with angular frequency $\omega$ and velocity potential

$$
\Phi(x,y,z,t) = \textrm{Re}( \phi(x,y,z)\textrm{e}^{-\textrm{i}\omega t} )
$$

is assumed.

In this expression, $\textrm{Re}( )$ denotes the real part of the complex expression, $i=\sqrt{-1}$, $t$ represents time, and $\phi(x,y,z)$ is the complex amplitude of the total velocity potential.

### Governing Equations and Boundary Conditions

Under the previous assumptions, the 3D Laplace equation governs the flow field:

$$
\left(
\frac{\partial^{2}}{\partial x^{2}}+
\frac{\partial^{2}}{\partial y^{2}}+
\frac{\partial^{2}}{\partial z^{2}}
\right) \phi=0,
$$

The linearity of the problem allows decomposition of the velocity potential as

$$
\phi = \phi^{I} + \phi^{S},
$$

where $\phi^{I}$ and $\phi^{S}$ are the incident and scattered velocity potentials, respectively.

The incident velocity potential is

$$
\phi^{I} = -\frac{\textrm{i}g A}{\omega}
\frac{\cosh k ( z + h )}{\cosh k h}
\textrm{e}^{ikx}
$$

where $A$ is the incident wave amplitude, $g$ is gravitational acceleration and $k$ is the wave number obtained from the dispersion relation

$$
\omega^{2}=gk\tanh(kh)
$$

### Free Surface Condition

At the still water level, the linearized free-surface condition is

$$
\frac{\partial \phi}{\partial z}- K\phi = 0,
\quad (x,y,z)\in S_{FS},
$$

where $K=\omega^{2}/g$.


### Seabed and Back Wall Boundary Conditions

The seabed and wave flume walls are assumed impermeable. Therefore,

Bottom boundary:

$$
\frac{\partial \phi}{\partial z}=0,
\quad (x,y,z)\in S_{BOTTOM}
$$

Back wall:

$$
\frac{\partial \phi}{\partial x}=0,
\quad (x,y,z)\in S_{FRONT}, S_{BACK}
$$

### Radiation Condition

At the far field, the Sommerfeld radiation condition is applied:

$$
\frac{\partial (\phi - \phi^{I})}{\partial x} + i k (\phi - \phi^{I}) = 0, \qquad x \to -\infty
$$

$$
\frac{\partial \phi}{\partial x} + i k \phi  = 0, \qquad x \to +\infty
$$

### Continuity at Domain Interfaces (Multi-domain)

Continuity of potential and flux is enforced at the interface of Domain 0 and Domain 1 ($S_{OUT,0}$ and $S_{IN,1}$):

$$
\phi_{0} = \phi_{1}, \qquad 
\vec{n_{0}}\cdot\nabla\phi_{0} = -\vec{n_{1}}\cdot\nabla\phi_{1}
$$

The gradient operator is

$$
\nabla =
\frac{\partial}{\partial x}\mathbf{i}+
\frac{\partial}{\partial y}\mathbf{j}+
\frac{\partial}{\partial z}\mathbf{k}
$$

The minus sign indicates that the outgoing flux from one subdomain equals the incoming flux in the adjacent domain because outward normals are opposite.


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
- numpy
- scipy
- matplotlib
- gmsh (Python API)

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


