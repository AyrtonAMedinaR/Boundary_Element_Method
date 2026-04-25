// ============================================================
// GMSH PROJECT
// Rectangular outlet domain located downstream of cylinder
// ============================================================

//+SetFactory("OpenCASCADE");
//+Geometry.OCCAutoFix = 0;


// ============================================================
// GLOBAL MESH SETTINGS
// ============================================================

lc  = 0.5;                       // Characteristic mesh length

Mesh.ElementOrder = 2;           // Quadratic elements
Mesh.SecondOrderLinear = 0;      // Curved second-order geometry


// ============================================================
// GEOMETRIC PARAMETERS
// ============================================================

// --- Water depth ---
hd         = 0.50;

// --- Cylinder / structure parameters ---
RADIUS     = 0.45;               // Cylinder radius
WIDTH      = 0.10*hd;            // Channel width (y-direction)

// --- Longitudinal spacing ---
DISTANCE_X  = 0.10*hd;           // Distance to the interface (unused here)

// --- Outlet domain position and length ---
DOMAIN_0 = 3.00*hd + 2.00*RADIUS + 2*DISTANCE_X;   // Start of outlet domain
LENGTH_0 = 3.00*hd;                                // Outlet length


// ============================================================
// TRANFINITE DISCRETIZATION PARAMETERS
// ============================================================

ELEM       = 11;      // Vertical discretization
ELEM_ARCS  = 11;      // Height discretization


// ============================================================
// DOMAIN CORNER POINTS
// Rectangular channel volume (2D surfaces)
// ============================================================

// --- Bottom plane (z = -hd)
Point(1) = {DOMAIN_0 + 0,        -WIDTH/2, -hd, lc};
Point(2) = {DOMAIN_0 + LENGTH_0, -WIDTH/2, -hd, lc};
Point(3) = {DOMAIN_0 + 0,        +WIDTH/2, -hd, lc};
Point(4) = {DOMAIN_0 + LENGTH_0, +WIDTH/2, -hd, lc};

// --- Free surface (z = 0)
Point(5) = {DOMAIN_0 + 0,        -WIDTH/2, 0, lc};
Point(6) = {DOMAIN_0 + LENGTH_0, -WIDTH/2, 0, lc};
Point(7) = {DOMAIN_0 + 0,        +WIDTH/2, 0, lc};
Point(8) = {DOMAIN_0 + LENGTH_0, +WIDTH/2, 0, lc};


// ============================================================
// DOMAIN EDGES
// ============================================================

// --- Spanwise edges (Y direction)
Line(1) = {5, 7};
Line(2) = {1, 3};
Line(3) = {2, 4};
Line(4) = {6, 8};

// --- Streamwise edges (X direction)
Line(5) = {6, 5};
Line(6) = {2, 1};
Line(7) = {4, 3};
Line(8) = {8, 7};

// --- Vertical edges (Z direction)
Line(9)  = {1, 5};
Line(10) = {3, 7};
Line(11) = {2, 6};
Line(12) = {4, 8};


// ============================================================
// SURFACE DEFINITIONS
// ============================================================

// --- Entrance plane (upstream face of this block)
Curve Loop(1) = {1, -10, -2, 9};
Plane Surface(1) = {1};

// --- Outlet plane (downstream boundary)
Curve Loop(2) = {3, 12, -4, -11};
Plane Surface(2) = {2};

// --- Bottom surface
Curve Loop(3) = {5, -9, -6, 11};
Plane Surface(3) = {3};

// --- Top surface (free surface)
Curve Loop(4) = {7, 10, -8, -12};
Plane Surface(4) = {4};

// --- Lateral vertical surface
Curve Loop(5) = {8, -1, -5, 4};
Plane Surface(5) = {5};

// (Opposite lateral face intentionally omitted)


// ============================================================
// STRUCTURED MESH
// ============================================================

Transfinite Surface {1:6};
Recombine Surface {1:6};

// --- Streamwise discretization (X)
Transfinite Line {5:8} = 3*(ELEM - 1) + 1 Using Progression 1;

// --- Spanwise discretization (Y)
Transfinite Line {1:4} = 2 Using Progression 1;

// --- Vertical discretization (Z)
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

Mesh 2;   // Generate 2D surface mesh