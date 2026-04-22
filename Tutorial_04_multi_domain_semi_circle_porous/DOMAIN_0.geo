// ============================================================
// GMSH PROJECT
// Created: Wed Dec 14 13:30:10 2022
// Kernel
// ============================================================

//+SetFactory("OpenCASCADE");
//+Geometry.OCCAutoFix = 0;


// ============================================================
// MESH GLOBAL SETTINGS
// ============================================================

lc  = 0.5;                       // Characteristic length (global mesh size)

Mesh.ElementOrder = 2;           // Quadratic elements
Mesh.SecondOrderLinear = 0;      // Use curved second-order elements


// ============================================================
// GEOMETRIC PARAMETERS (CHANNEL / DOMAIN)
// ============================================================

// --- Vertical dimension ---
hd         = 0.50;               // Water depth (domain height)

// --- Cylinder / structure parameters ---
RADIUS     = 0.45;               // Cylinder radius
WIDTH      = 0.10*hd;            // Channel width (y-direction)

// --- Longitudinal spacing ---
DISTANCE_X  = 0.10*hd;           // Distance to the interface (unused here)

// --- Domain limits in X ---
DOMAIN_0 = 0;                    // Inlet x-position
LENGTH_0 = 3.00*hd;              // Domain length


// ============================================================
// TRANFINITE DISCRETIZATION PARAMETERS
// ============================================================

ELEM = 11;                       // Base discretization parameter
ELEM_ARCS = 11;                  // Vertical discretization


// ============================================================
// DOMAIN CORNER POINTS
// Rectangular channel: (x, y, z)
// z = -hd  -> bottom
// z =  0   -> free surface
// ============================================================

// Bottom plane (z = -hd)
Point(1) =  {DOMAIN_0 + 0,        -WIDTH/2, -hd,  lc};
Point(2) =  {DOMAIN_0 + LENGTH_0, -WIDTH/2, -hd,  lc};
Point(3) =  {DOMAIN_0 + 0,        +WIDTH/2, -hd,  lc};
Point(4) =  {DOMAIN_0 + LENGTH_0, +WIDTH/2, -hd,  lc};

// Free surface (z = 0)
Point(5) =  {DOMAIN_0 + 0,        -WIDTH/2,  0,  lc};
Point(6) =  {DOMAIN_0 + LENGTH_0, -WIDTH/2,  0,  lc};
Point(7) =  {DOMAIN_0 + 0,        +WIDTH/2,  0,  lc};
Point(8) =  {DOMAIN_0 + LENGTH_0, +WIDTH/2,  0,  lc};


// ============================================================
// EDGE LINES OF THE CHANNEL BOX
// ============================================================

// Top edges (free surface rectangle)
Line(1) = {5, 7};
Line(4) = {6, 8};

// Bottom edges
Line(2) = {1, 3};
Line(3) = {2, 4};

// Longitudinal edges (x-direction)
Line(5) = {6, 5};
Line(6) = {2, 1};
Line(7) = {4, 3};
Line(8) = {8, 7};

// Vertical edges (z-direction)
Line(9)  = {1, 5};
Line(10) = {3, 7};
Line(11) = {2, 6};
Line(12) = {4, 8};


// ============================================================
// SURFACE DEFINITION (CHANNEL FACES)
// ============================================================

// Inlet surface (x = DOMAIN_0)
Curve Loop(1) = {1, -10, -2, 9};
Plane Surface(1) = {1};

// Outlet surface (x = DOMAIN_0 + LENGTH_0)
Curve Loop(2) = {3, 12, -4, -11};
Plane Surface(2) = {2};

// Front wall (y = -WIDTH/2)
Curve Loop(3) = {5, -9, -6, 11};
Plane Surface(3) = {3};

// Back wall (y = +WIDTH/2)
Curve Loop(4) = {7, 10, -8, -12};
Plane Surface(4) = {4};

// Free surface (z = 0)
Curve Loop(5) = {8, -1, -5, 4};
Plane Surface(5) = {5};

// Bottom surface (z = -hd)  (currently commented)
//+Curve Loop(6) = {6, 2, -7, -3};
//+Plane Surface(6) = {6};


// ============================================================
// TRANSFINITE STRUCTURED MESH
// ============================================================

Transfinite Surface {1:6};       // Structured surfaces
Recombine Surface {1:6};         // Quad elements

// --- Discretization along X direction ---
Transfinite Line {5:8} = 3*(ELEM - 1) + 1 Using Progression 1;

// --- Discretization along Y direction ---
Transfinite Line {1:4} = 2 Using Progression 1;

// --- Discretization along Z (depth) ---
Transfinite Line {9:12} = ELEM_ARCS Using Progression 1;

Mesh.Smoothing = 1;


// ============================================================
// PHYSICAL GROUPS (BOUNDARY CONDITIONS)
// ============================================================

Physical Surface("ENTRANCE", 1) = {1};
Physical Surface("OUTLET",2) = {2};

Physical Surface("POROUS_BOUNDARY",3) = { };
Physical Surface("FREE_SURFACE",4) = {5};
Physical Surface("BOTTOM",5) = { };
Physical Surface("IMPERMEABLE",6) = {};

Physical Surface("FRONT",7) = {3};
Physical Surface("BACK",8) = {4};


// ============================================================
// MESH GENERATION
// ============================================================

Mesh 2;   // Generate 2D mesh of surfaces