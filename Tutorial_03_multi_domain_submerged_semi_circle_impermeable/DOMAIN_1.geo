// ============================================================
// Gmsh geometry file
// Domain with half-cylinder obstacle
// ============================================================

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
hd         = 0.50;              // Water depth (Z)
WIDTH      = 0.25 * hd;         // Width (Y)
TRANSLATION   = 3.00 * hd;         // X start (continuation of previous domain)

// Cylinder parameters
RADIUS     = 0.45;
DX = 0.10 * hd;

// Mesh resolution
ELEM       = 7;
ELEM_ARCS  = 6;

lc = hd/10;

// ------------------------------------------------------------
// POINTS
// ------------------------------------------------------------

// Bottom plane (Z = -hd)
Point(1) = {TRANSLATION, -WIDTH/2, -hd, lc};
Point(2) = {TRANSLATION + 2*DX + 2*RADIUS, -WIDTH/2, -hd, lc};
Point(3) = {TRANSLATION,  WIDTH/2, -hd, lc};
Point(4) = {TRANSLATION + 2*DX + 2*RADIUS,  WIDTH/2, -hd, lc};

// Free surface (Z = 0)
Point(5) = {TRANSLATION, -WIDTH/2, 0, lc};
Point(6) = {TRANSLATION + 2*DX + 2*RADIUS, -WIDTH/2, 0, lc};
Point(7) = {TRANSLATION,  WIDTH/2, 0, lc};
Point(8) = {TRANSLATION + 2*DX + 2*RADIUS,  WIDTH/2, 0, lc};

// Half-cylinder (bottom edge, Y = -W/2)
Point(9)  = {TRANSLATION + DX,              -WIDTH/2, -hd, lc};
Point(10) = {TRANSLATION + DX + RADIUS,     -WIDTH/2, -hd, lc};
Point(11) = {TRANSLATION + DX + 2*RADIUS,   -WIDTH/2, -hd, lc};

Point(12) = {TRANSLATION + DX + RADIUS - RADIUS*Cos(45*Pi/180),
             -WIDTH/2,
             -hd + RADIUS*Cos(45*Pi/180), lc};

Point(13) = {TRANSLATION + DX + RADIUS + RADIUS*Cos(45*Pi/180),
             -WIDTH/2,
             -hd + RADIUS*Cos(45*Pi/180), lc};

// Half-cylinder (top edge, Y = +W/2)
Point(14) = {TRANSLATION + DX,              WIDTH/2, -hd, lc};
Point(15) = {TRANSLATION + DX + RADIUS,     WIDTH/2, -hd, lc};
Point(16) = {TRANSLATION + DX + 2*RADIUS,   WIDTH/2, -hd, lc};

Point(17) = {TRANSLATION + DX + RADIUS - RADIUS*Cos(45*Pi/180),
             WIDTH/2,
             -hd + RADIUS*Cos(45*Pi/180), lc};

Point(18) = {TRANSLATION + DX + RADIUS + RADIUS*Cos(45*Pi/180),
             WIDTH/2,
             -hd + RADIUS*Cos(45*Pi/180), lc};

// ------------------------------------------------------------
// LINES
// ------------------------------------------------------------

// Domain edges
Line(1) = {5, 7};
Line(2) = {1, 3};
Line(7) = {2, 4};
Line(8) = {6, 8};

Line(9)  = {1, 5};
Line(10) = {3, 7};
Line(11) = {2, 6};
Line(12) = {4, 8};

Line(13) = {6, 5};
Line(14) = {8, 7};

// Connections around cylinder
Line(3) = {9, 14};
Line(4) = {12, 17};
Line(5) = {13, 18};
Line(6) = {11, 16};

Line(15) = {2, 11};
Line(16) = {9, 1};
Line(17) = {4, 16};
Line(18) = {14, 3};

Line(19) = {12, 5};
Line(20) = {13, 6};
Line(21) = {17, 7};
Line(22) = {18, 8};

// ------------------------------------------------------------
// HALF-CYLINDER ARCS
// ------------------------------------------------------------

Circle(23) = {9, 10, 12};
Circle(24) = {12, 10, 13};
Circle(25) = {13, 10, 11};

Circle(26) = {16, 15, 18};
Circle(27) = {18, 15, 17};
Circle(28) = {17, 15, 14};

// ------------------------------------------------------------
// SURFACES
// ------------------------------------------------------------

// Domain boundaries
Curve Loop(1) = {1, -10, -2, 9};
Plane Surface(1) = {1};   // Entrance

Curve Loop(2) = {7, 12, -8, -11};
Plane Surface(2) = {2};   // Outlet

Curve Loop(3) = {14, -1, -13, 8};
Plane Surface(3) = {3};   // Free surface

// Cylinder side patches (front Y-)
Curve Loop(6) = {23, 19, -9, -16};
Plane Surface(6) = {6};

Curve Loop(7) = {24, 20, 13, -19};
Plane Surface(7) = {7};

Curve Loop(8) = {25, -15, 11, -20};
Plane Surface(8) = {8};

// Cylinder side patches (back Y+)
Curve Loop(9) = {26, 22, -12, 17};
Plane Surface(9) = {9};

Curve Loop(10) = {27, 21, -14, -22};
Plane Surface(10) = {10};

Curve Loop(11) = {28, 18, 10, -21};
Plane Surface(11) = {11};

// Cylinder curved surface (split into 3)
Curve Loop(12) = {3, -28, -4, -23};
Surface(12) = {12};

Curve Loop(14) = {4, -27, -5, -24};
Surface(13) = {14};

Curve Loop(16) = {5, -26, -6, -25};
Surface(14) = {16};

// ------------------------------------------------------------
// STRUCTURED MESH
// ------------------------------------------------------------

Transfinite Surface {1:14};
Recombine Surface {1:14};

// Cylinder arcs
Transfinite Line {24, 27} = ELEM Using Progression 1;
Transfinite Line {23, 25, 26, 28} = ELEM_ARCS Using Progression 1;

// Near-cylinder refinement
Transfinite Line {15:18} = 5 Using Progression 1;
Transfinite Line {19:22} = 5 Using Progression 1;

// Domain resolution
Transfinite Line {13:14} = ELEM Using Progression 1;
Transfinite Line {1:8}   = 3 Using Progression 1;
Transfinite Line {9:12}  = ELEM_ARCS Using Progression 1;

// ------------------------------------------------------------
// PHYSICAL GROUPS
// ------------------------------------------------------------

Physical Surface("ENTRANCE", 1) = {1};

Physical Surface("OUTLET",2) = {2};

Physical Surface("POROUS_BOUNDARY",3) = {};

Physical Surface("FREE_SURFACE",4) = {3};

Physical Surface("BOTTOM",5) = { };

Physical Surface("IMPERMEABLE",6) = {12:14};

Physical Surface("FRONT",7) = {6:8};

Physical Surface("BACK",8) = {9:11};

// ------------------------------------------------------------
// MESH GENERATION
// ------------------------------------------------------------
Mesh 2;
