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

// --- Longitudinal spacing ---
DISTANCE_X  = 0.10*hd;           // Gap between cylinder and boundaries

// --- Domain location and length ---
DOMAIN_0 = 3.00*hd;              // Start of this block in x
LENGTH_0 = 3.00*hd;              // Length of this block


// ============================================================
// TRANFINITE DISCRETIZATION PARAMETERS
// ============================================================

ELEM = 11;                       // Main discretization parameter
ELEM_ARCS = 11;                  // Arc / vertical discretization


// ============================================================
// DOMAIN CORNER POINTS
// Rectangular section containing the half-cylinder
// ============================================================

// --------------------
// Bottom plane (z = -hd)
// --------------------
Point(1) =  {DOMAIN_0 + 0, -WIDTH/2, -hd,  lc};
Point(2) =  {DOMAIN_0 + 2*DISTANCE_X + 2*RADIUS, -WIDTH/2, -hd,  lc};
Point(3) =  {DOMAIN_0 + 0, WIDTH/2, -hd,  lc};
Point(4) =  {DOMAIN_0 + 2*DISTANCE_X + 2*RADIUS, WIDTH/2, -hd,  lc};

// --------------------
// Free surface (z = 0)
// --------------------
Point(5) =  {DOMAIN_0 + 0, -WIDTH/2, 0,  lc};
Point(6) =  {DOMAIN_0 + 2*DISTANCE_X + 2*RADIUS, -WIDTH/2, 0,  lc};
Point(7) =  {DOMAIN_0 + 0, WIDTH/2, 0,  lc};
Point(8) =  {DOMAIN_0 + 2*DISTANCE_X + 2*RADIUS, WIDTH/2, 0,  lc};


// ============================================================
// HALF-CYLINDER GEOMETRY (BOTTOM ATTACHED)
// Points defining two semicircular arcs (front/back)
// ============================================================

// --- Front side (y = -WIDTH/2)
Point(9)  = {DOMAIN_0 + DISTANCE_X, -WIDTH/2, -hd,  lc};
Point(10) = {DOMAIN_0 + DISTANCE_X + RADIUS, -WIDTH/2, -hd,  lc};
Point(11) = {DOMAIN_0 + DISTANCE_X + 2*RADIUS, -WIDTH/2, -hd,  lc};
Point(12) = {DOMAIN_0 + DISTANCE_X + RADIUS - RADIUS*Cos(45*Pi/180), -WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};
Point(13) = {DOMAIN_0 + DISTANCE_X + RADIUS + RADIUS*Cos(45*Pi/180), -WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};

// --- Back side (y = +WIDTH/2)
Point(14) = {DOMAIN_0 + DISTANCE_X, +WIDTH/2, -hd,  lc};
Point(15) = {DOMAIN_0 + DISTANCE_X + RADIUS, +WIDTH/2, -hd,  lc};
Point(16) = {DOMAIN_0 + DISTANCE_X + 2*RADIUS, +WIDTH/2, -hd,  lc};
Point(17) = {DOMAIN_0 + DISTANCE_X + RADIUS - RADIUS*Cos(45*Pi/180), +WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};
Point(18) = {DOMAIN_0 + DISTANCE_X + RADIUS + RADIUS*Cos(45*Pi/180), +WIDTH/2, -hd + RADIUS*Cos(45*Pi/180),  lc};


// ============================================================
// STRAIGHT EDGES (LINES)
// Channel box + cylinder connectors
// ============================================================

Line(1) = {5, 7};
Line(2) = {1, 3};
Line(3) = {9, 14};
Line(4) = {12, 17};
Line(5) = {13, 18};
Line(6) = {11, 16};
Line(7) = {2, 4};
Line(8) = {6, 8};

Line(9)  = {1, 5};
Line(10) = {3, 7};
Line(11) = {2, 6};
Line(12) = {4, 8};

Line(13) = {6, 5};
Line(14) = {8, 7};

Line(15) = {2, 11};
Line(16) = {9, 1};
Line(17) = {4, 16};
Line(18) = {14, 3};

Line(19) = {12, 5};
Line(20) = {13, 6};
Line(21) = {17, 7};
Line(22) = {18, 8};


// ============================================================
// CYLINDER ARCS (SEMICIRCLES FRONT/BACK)
// ============================================================

Circle(23) = {9, 10, 12};
Circle(24) = {12, 10, 13};
Circle(25) = {13, 10, 11};

Circle(26) = {16, 15, 18};
Circle(27) = {18, 15, 17};
Circle(28) = {17, 15, 14};


// ============================================================
// SURFACE DEFINITIONS
// ============================================================

// Inlet
Curve Loop(1) = {1, -10, -2, 9};
Plane Surface(1) = {1};

// Outlet
Curve Loop(2) = {7, 12, -8, -11};
Plane Surface(2) = {2};

// Free surface block
Curve Loop(3) = {14, -1, -13, 8};
Plane Surface(3) = {3};

// Half-cylinder front patches
Curve Loop(6) = {23, 19, -9, -16};
Plane Surface(6) = {6};

Curve Loop(7) = {24, 20, 13, -19};
Plane Surface(7) = {7};

Curve Loop(8) = {25, -15, 11, -20};
Plane Surface(8) = {8};

// Half-cylinder back patches
Curve Loop(9)  = {26, 22, -12, 17};
Plane Surface(9)  = {9};

Curve Loop(10) = {27, 21, -14, -22};
Plane Surface(10) = {10};

Curve Loop(11) = {28, 18, 10, -21};
Plane Surface(11) = {11};

// Side porous faces (extruded surfaces)
Curve Loop(12) = {3, -28, -4, -23};
Surface(12) = {12};

Curve Loop(14) = {4, -27, -5, -24};
Surface(13) = {14};

Curve Loop(16) = {5, -26, -6, -25};
Surface(14) = {16};


// ============================================================
// STRUCTURED MESH (TRANSFINITE + RECOMBINE)
// ============================================================

Transfinite Surface {1:14};
Recombine Surface {1:14};


// --- Cylinder arcs ---
Transfinite Line {24,27} = ELEM Using Progression 1;
Transfinite Line {23,25,26,28} = ELEM_ARCS Using Progression 1;

// --- Squares around cylinder ---
Transfinite Line {15:18} = 5 Using Progression 1;
Transfinite Line {19:22} = 5 Using Progression 1;

// --- Streamwise direction ---
Transfinite Line {13:14} = ELEM Using Progression 1;

// --- Spanwise direction ---
Transfinite Line {1:8} = 2 Using Progression 1;

// --- Vertical direction ---
Transfinite Line {9:12} = ELEM_ARCS Using Progression 1;

Mesh.Smoothing = 1;


// ============================================================
// PHYSICAL GROUPS (BOUNDARY CONDITIONS)
// ============================================================

Physical Surface("ENTRANCE", 1) = {1};
Physical Surface("OUTLET",2) = {2};

Physical Surface("POROUS_BOUNDARY",3) = {12:14};
Physical Surface("FREE_SURFACE",4) = {3};
Physical Surface("BOTTOM",5) = { };
Physical Surface("IMPERMEABLE",6) = {};

Physical Surface("FRONT",7) = {6:8};
Physical Surface("BACK",8) = {9:11};


// ============================================================
// MESH GENERATION
// ============================================================

Mesh 2;   // Generate 2D surface mesh