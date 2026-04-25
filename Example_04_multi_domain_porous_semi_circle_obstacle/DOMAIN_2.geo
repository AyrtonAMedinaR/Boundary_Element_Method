// ============================================================
// GMSH PROJECT
// Created: Wed Dec 14 13:30:10 2022
// ============================================================

//+SetFactory("OpenCASCADE");
//+Geometry.OCCAutoFix = 0;


// ============================================================
// GLOBAL MESH SETTINGS
// ============================================================

lc  = 0.5;                       // Global characteristic length

Mesh.ElementOrder = 2;           // Quadratic elements
Mesh.SecondOrderLinear = 0;      // Curved second-order geometry


// ============================================================
// GEOMETRIC PARAMETERS
// ============================================================

// --- Vertical scale ---
hd         = 0.50;               // Water depth

// --- Cylinder parameters ---
RADIUS     = 0.45;               // Cylinder radius
WIDTH      = 0.10*hd;            // Channel width (y-direction)

// --- Streamwise offset of the cylinder block ---
DX  = 0.10*hd;                   // Distance from inlet block
LX = 3.00*hd;                    // Location of this geometry block


// ============================================================
// TRANFINITE DISCRETIZATION PARAMETERS
// ============================================================

ELEM = 11;                       // Main discretization parameter
ELEM_ARCS = 11;                  // Arc / vertical discretization


// ============================================================
// HALF-CYLINDER GEOMETRY (BOTTOM ATTACHED)
// Two semicircular faces + interior square refinement
// ============================================================

// ------------------------------------------------------------
// MAIN ARC POINTS (front and back faces of cylinder)
// ------------------------------------------------------------

// --- Front arc (y = -WIDTH/2)
Point(1) =  {LX + DX, -WIDTH/2, -hd,  lc};
Point(2) =  {LX + DX + RADIUS, -WIDTH/2, -hd,  lc};
Point(3) =  {LX + DX + 2*RADIUS, -WIDTH/2, -hd,  lc};
Point(4) =  {LX + DX + RADIUS - RADIUS*Cos(45*Pi/180), -WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};
Point(5) =  {LX + DX + RADIUS + RADIUS*Cos(45*Pi/180), -WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};

// --- Back arc (y = +WIDTH/2)
Point(6)  = {LX + DX, +WIDTH/2, -hd,  lc};
Point(7)  = {LX + DX + RADIUS, +WIDTH/2, -hd,  lc};
Point(8)  = {LX + DX + 2*RADIUS, +WIDTH/2, -hd,  lc};
Point(9)  = {LX + DX + RADIUS - RADIUS*Cos(45*Pi/180), +WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};
Point(10) = {LX + DX + RADIUS + RADIUS*Cos(45*Pi/180), +WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};


// ------------------------------------------------------------
// INTERNAL REFINEMENT POINTS (square patch inside cylinder)
// ------------------------------------------------------------

// --- Front interior square
Point(11) = {LX + DX + RADIUS/2,   -WIDTH/2, -hd,  lc};
Point(12) = {LX + DX + 3*RADIUS/2, -WIDTH/2, -hd,  lc};
Point(13) = {LX + DX + RADIUS/2,   -WIDTH/2, -hd + RADIUS/2,  lc};
Point(14) = {LX + DX + 3*RADIUS/2, -WIDTH/2, -hd + RADIUS/2,  lc};

// --- Back interior square
Point(15) = {LX + DX + RADIUS/2,   WIDTH/2, -hd,  lc};
Point(16) = {LX + DX + 3*RADIUS/2, WIDTH/2, -hd,  lc};
Point(17) = {LX + DX + RADIUS/2,   WIDTH/2, -hd + RADIUS/2,  lc};
Point(18) = {LX + DX + 3*RADIUS/2, WIDTH/2, -hd + RADIUS/2,  lc};


// ============================================================
// CYLINDER ARCS (SEMICIRCLES)
// ============================================================

Circle(1) = {1, 2, 4};
Circle(2) = {4, 2, 5};
Circle(3) = {5, 2, 3};

Circle(4) = {6, 7, 9};
Circle(5) = {9, 7, 10};
Circle(6) = {10, 7, 8};


// ============================================================
// CONNECTING LINES
// ============================================================

// --- Spanwise connections (front ↔ back)
Line(7)  = {1, 6};
Line(8)  = {4, 9};
Line(9)  = {5, 10};
Line(10) = {3, 8};


// ------------------------------------------------------------
// FRONT PATCH (refinement square)
// ------------------------------------------------------------

// Diagonals
Line(11) = {1, 11};
Line(12) = {4, 13};
Line(13) = {5, 14};
Line(14) = {3, 12};

// Horizontal edges
Line(15) = {11, 12};
Line(16) = {13, 14};

// Vertical edges
Line(17) = {13, 11};
Line(18) = {14, 12};


// ------------------------------------------------------------
// BACK PATCH (refinement square)
// ------------------------------------------------------------

// Diagonals
Line(19) = {6, 15};
Line(20) = {9, 17};
Line(21) = {10, 18};
Line(22) = {8, 16};

// Horizontal edges
Line(23) = {15, 16};
Line(24) = {17, 18};

// Vertical edges
Line(25) = {17, 15};
Line(26) = {18, 16};


// ============================================================
// SURFACE DEFINITIONS
// ============================================================

// --- Curved cylindrical surfaces (porous boundary)
Curve Loop(1) = {8, -4, -7, 1};
Curve Loop(2) = {9, -5, -8, 2};
Curve Loop(3) = {10, -6, -9, 3};

Surface(1) = {1};
Surface(2) = {2};
Surface(3) = {3};


// --- Front faces
Curve Loop(4) = {11, -17, -12, -1};
Curve Loop(5) = {16, -13, -2, 12};
Curve Loop(6) = {13, 18, -14, -3};
Curve Loop(7) = {17, 15, -18, -16};

Plane Surface(4) = {4};
Plane Surface(5) = {5};
Plane Surface(6) = {6};
Plane Surface(7) = {7};


// --- Back faces
Curve Loop(8)  = {20, 25, -19, 4};
Curve Loop(9)  = {21, -24, -20, 5};
Curve Loop(10) = {22, -26, -21, 6};
Curve Loop(11) = {26, -23, -25, 24};

Plane Surface(8)  = {8};
Plane Surface(9)  = {9};
Plane Surface(10) = {10};
Plane Surface(11) = {11};


// ============================================================
// STRUCTURED MESH (TRANSFINITE + RECOMBINE)
// ============================================================

Transfinite Surface {1:11};
Recombine Surface {1:11};

// --- Arc discretization
Transfinite Line {2,5} = ELEM Using Progression 1;
Transfinite Line {1,3,4,6} = ELEM_ARCS Using Progression 1;

// --- Interior square discretization
Transfinite Line {15:16,23:24} = ELEM Using Progression 1;
Transfinite Line {17:18,25:26} = ELEM_ARCS Using Progression 1;

// --- Radial connectors
Transfinite Line {11:14,19:22} = 5 Using Progression 1;

// --- Spanwise direction
Transfinite Line {7:10} = 2 Using Progression 1;


// ============================================================
// PHYSICAL GROUPS (BOUNDARY CONDITIONS)
// ============================================================

Physical Surface("ENTRANCE", 1) = { };
Physical Surface("OUTLET",2) = { };

Physical Surface("POROUS_BOUNDARY",3) = {1:3};

Physical Surface("FREE_SURFACE",4) = { };
Physical Surface("BOTTOM",5) = { };
Physical Surface("IMPERMEABLE",6) = {};

Physical Surface("FRONT",7) = {4:7};
Physical Surface("BACK",8) = {8:11};


// ============================================================
// MESH GENERATION
// ============================================================

Mesh 2;   // Generate 2D surface mesh