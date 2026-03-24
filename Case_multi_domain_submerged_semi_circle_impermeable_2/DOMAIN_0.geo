// ============================================================
// Gmsh geometry file
// Rectangular 3D domain (extruded in Z)
// ============================================================

SetFactory("OpenCASCADE");

// ------------------------------------------------------------
// MESH SETTINGS
// ------------------------------------------------------------

Mesh.ElementOrder = 2;
Mesh.SecondOrderLinear = 0;
Mesh.Smoothing = 1;

// ------------------------------------------------------------
// GEOMETRY PARAMETERS
// ------------------------------------------------------------

// Domain dimensions
hd        = 0.40;           // Water depth (Z direction)
LENGTH    = 3.00 * hd;      // Length in X
WIDTH     = 0.40 * hd;     // Width in Y

// Mesh resolution parameters
ELEM       = 7;
ELEM_ARCS  = 6;

lc = hd/10;

// ------------------------------------------------------------
// POINTS (Bottom Z = -hd, Top Z = 0)
// ------------------------------------------------------------

// Bottom plane
Point(1) = {0,        -WIDTH/2, -hd, lc};
Point(2) = {LENGTH,   -WIDTH/2, -hd, lc};
Point(3) = {0,         WIDTH/2, -hd, lc};
Point(4) = {LENGTH,    WIDTH/2, -hd, lc};

// Free surface (top)
Point(5) = {0,        -WIDTH/2, 0, lc};
Point(6) = {LENGTH,   -WIDTH/2, 0, lc};
Point(7) = {0,         WIDTH/2, 0, lc};
Point(8) = {LENGTH,    WIDTH/2, 0, lc};

// ------------------------------------------------------------
// LINES (Edges of the box)
// ------------------------------------------------------------

// Y-direction edges
Line(1) = {5, 7};
Line(2) = {1, 3};
Line(3) = {2, 4};
Line(4) = {6, 8};

// X-direction edges
Line(5) = {6, 5};
Line(6) = {2, 1};
Line(7) = {4, 3};
Line(8) = {8, 7};

// Vertical edges (Z-direction)
Line(9)  = {1, 5};
Line(10) = {3, 7};
Line(11) = {2, 6};
Line(12) = {4, 8};

// ------------------------------------------------------------
// SURFACES
// ------------------------------------------------------------

// Entrance (X = 0)
Curve Loop(1) = {1, -10, -2, 9};
Plane Surface(1) = {1};

// Outlet (X = L)
Curve Loop(2) = {3, 12, -4, -11};
Plane Surface(2) = {2};

// Front (Y = -W/2)
Curve Loop(3) = {5, -9, -6, 11};
Plane Surface(3) = {3};

// Back (Y = +W/2)
Curve Loop(4) = {7, 10, -8, -12};
Plane Surface(4) = {4};

// Free surface (Z = 0)
Curve Loop(5) = {8, -1, -5, 4};
Plane Surface(5) = {5};

// ------------------------------------------------------------
// STRUCTURED MESH (TRANSFINITE)
// ------------------------------------------------------------

// Apply structured meshing
Transfinite Surface {1:5};
Recombine Surface {1:5};

// X-direction resolution
Transfinite Line {5:8} = 3*(ELEM - 1) + 1 Using Progression 1;

// Y-direction resolution
Transfinite Line {1:4} = 3 Using Progression 1;

// Z-direction resolution
Transfinite Line {9:12} = ELEM_ARCS Using Progression 1;

// ------------------------------------------------------------
// PHYSICAL GROUPS (BOUNDARY CONDITIONS)
// ------------------------------------------------------------

Physical Surface("ENTRANCE", 1) = {1};

Physical Surface("OUTLET",2) = {2};

Physical Surface("POROUS_BOUNDARY",3) = { };

Physical Surface("FREE_SURFACE",4) = {5};

Physical Surface("BOTTOM",5) = { };

Physical Surface("IMPERMEABLE",6) = {};

Physical Surface("FRONT",7) = {3};

Physical Surface("BACK",8) = {4};


// ------------------------------------------------------------
// MESH GENERATION
// ------------------------------------------------------------
Mesh 2;

