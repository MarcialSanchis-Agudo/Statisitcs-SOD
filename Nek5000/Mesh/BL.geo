// Obstacle
// - Geometry
H = 3.0;
W = 4.0;
W05 = W/2;
Lin = 10.0;
Lout = 7.5;

Lx = 0.37035082916161244;
Ly = 0.37035082916161244;
Ly_in = 0.37035082916161244;
Ly_out = 0.37035082916161244;
Ly_left = 0.37035082916161244;
Ly_right = 0.37035082916161244;
Ly_in_left = 0.37035082916161244;
Ly_in_right = 0.37035082916161244;
Ly_out_left = 0.37035082916161244;
Ly_out_right = 0.37035082916161244;

Lz = 0.37035082916161244;

Ne_wx = 17;
b_wx = 4.52;
Ne_wz = 6;
b_wz = 4.52;

Ne_xz = 19;
r_xz = 1.13;

Ne_y1 = 40;
r_y1 = 0.8071067811865486;
b_y1 = 2.8284271247461943;

Ne_x2_in = 43;
Ne_x2_out = 32;
Ne_y2 = 19;
r_y2 = 1;

Ne_z2 = 24;
r_z2 = 1.0;


// Exterior
// - Bottom Points:
Point(9) = {-Lin + 0,  0 + 0,  -W05 + -0.0};
Point(12) = {-Lin + 0,  0 + 0,  W05 + -0.0};


Point(18) = {Lout + 0,  0 + 0,  -W05 + -0.0};
Point(15) = {Lout + 0,  0 + 0,  W05 + -0.0};


// - Top Points:
Point(109) = {-Lin + 0,  H + 0,  -W05 + -0.0};
Point(112) = {-Lin + 0,  H + 0,  W05 + -0.0};

Point(118) = {Lout + 0,  H + 0,  -W05 + -0.0};
Point(115) = {Lout + 0,  H + 0,  W05 + -0.0};

// - Bottom z-Lines:
Line(13) = {9, 12};
Transfinite Curve{13} = Ne_z2 Using Progression r_z2;

Line(20) = {18, 15};
Transfinite Curve{20} = Ne_z2 Using Progression r_z2;

// - Top z-Lines:

Line(113) = {109, 112};
Transfinite Curve{113} = Ne_z2 Using Progression r_z2;

Line(120) = {118, 115};
Transfinite Curve{120} = Ne_z2 Using Progression r_z2;

// - Bottom x-lines:
Line(23) = {9, 18};
Transfinite Curve{23} = Ne_x2_in Using Progression 1;


Line(26) = {12, 15};
Transfinite Curve{26} = Ne_x2_in Using Progression 1;

// - Top x-lines:

Line(123) = {109, 118};
Transfinite Curve{123} = Ne_x2_in Using Progression 1;

Line(126) = {112, 115};
Transfinite Curve{126} = Ne_x2_in Using Progression 1;

// - Bottom Vertical lines
Line(209) = {9, 109};
Transfinite Curve{-209} = Ne_y1 Using Progression r_y1;

Line(212) = {12, 112};
Transfinite Curve{-212} = Ne_y1 Using Progression r_y1;

Line(215) = {15, 115};
Transfinite Curve{-215} = Ne_y1 Using Progression r_y1;

Line(218) = {18, 118};
Transfinite Curve{-218} = Ne_y1 Using Progression r_y1;


// - Bottom xz-plane
Line Loop(12) = {13, 26, -20, -23};
Plane Surface(12) = {12};

// - Bottom xz-plane
Line Loop(31) = {113, 126, -120, -123};
Plane Surface(31) = {31};

// - Bottom Vertical Surfaces:
Line Loop(113) = {-13, 209, 113, -212};
Plane Surface(113) = {113};

Line Loop(120) = {-20, 218, 120, -215};
Plane Surface(120) = {120};


Line Loop(123) = {23, 218, -123, -209};
Plane Surface(123) = {123};


Line Loop(126) = {26, 215, -126, -212};
Plane Surface(126) = {126};




// - Solids
// - - Bottom:
Surface Loop(12) = {12, 31, 113, 120, 123, 126};
Volume(12) = {12};




// Boundaries:
// Wall
// - Bottom: 2, 3, 4, 5, 12, 13, 14, 15, 16, 17, 18, 19
// - Top: 31, 32, 33, 34, 35, 36, 37, 38, 39
// - Left: 123, 127, 129, 173, 177, 179
// - Right: 126, 128, 132, 176, 178, 182
// - Obstacle: 101, 102, 103, 104, 6
// Inlet: 113, 114, 115, 163, 164, 165
// Outlet: 120, 121, 122, 170, 171, 172


Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";
Mesh 3;
SetOrder 2;
Physical Surface("inlet", 1) = {113};
Physical Surface("outlet", 2) = {120};
Physical Surface("bottom", 3) = {12};
Physical Surface("top", 4) = {31};
Physical Surface("left", 5) = {123};
Physical Surface("right", 6) = {126};
Physical Volume("fluid", 8) = {12};






