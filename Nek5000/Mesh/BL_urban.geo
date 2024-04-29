// Obstacle
// - Geometry
H = 3.0;
W = 4.0;
W05 = W/2;
Lin = 10.0;
Lout = 6;

Lx = 0.37035082916161244;
Ly = 1.3703508291616124;
Ly_in = 1.3703508291616124;
Ly_out = 1.3703508291616124;
Ly_left = 1.3703508291616124;
Ly_right = 1.3703508291616124;
Ly_in_left = 1.3703508291616124;
Ly_in_right = 1.3703508291616124;
Ly_out_left = 1.3703508291616124;
Ly_out_right = 1.3703508291616124;

Lz = 0.37035082916161244;

Ne_wx = 17;
b_wx = 4.52;
Ne_wz = 6;
b_wz = 4.52;

Ne_xz = 19;
r_xz = 1.13;

Ne_y1 = 31;
r_y1 = 1.104305130887847;
b_y1 = 4.417220523551388;

Ne_x2_in = 73;
Ne_x2_out = 32;
Ne_y2 = 12;
r_y2 = 1;

Ne_z2 = 27;
r_z2 = 1.0;


// Exterior
// - Bottom Points:
Point(9) = {-Lin + 0,  0 + 0,  -W05 + -0.0};
Point(12) = {-Lin + 0,  0 + 0,  W05 + -0.0};


Point(18) = {Lout + 0,  0 + 0,  -W05 + -0.0};
Point(15) = {Lout + 0,  0 + 0,  W05 + -0.0};

// - Middle points

Point(59) = {-Lin + 0,  Ly + 0,  -W05 + -0.0};
Point(62) = {-Lin + 0,  Ly + 0,  W05 + -0.0};


Point(68) = {Lout + 0,  Ly + 0,  -W05 + -0.0};
Point(65) = {Lout + 0,  Ly + 0,  W05 + -0.0};


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

// - Middle z-Lines:
Line(63) = {59, 62};
Transfinite Curve{63} = Ne_z2 Using Progression r_z2;

Line(70) = {68, 65};
Transfinite Curve{70} = Ne_z2 Using Progression r_z2;

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

// - Middle x-lines:
Line(73) = {59, 68};
Transfinite Curve{73} = Ne_x2_in Using Progression 1;


Line(76) = {62, 65};
Transfinite Curve{76} = Ne_x2_in Using Progression 1;

// - Top x-lines:

Line(123) = {109, 118};
Transfinite Curve{123} = Ne_x2_in Using Progression 1;

Line(126) = {112, 115};
Transfinite Curve{126} = Ne_x2_in Using Progression 1;

// - Bottom Vertical lines
Line(159) = {9, 59};
Transfinite Curve{-159} = Ne_y1 Using Progression 1/r_y1;

Line(162) = {12, 62};
Transfinite Curve{-162} = Ne_y1 Using Progression 1/r_y1;

Line(165) = {15, 65};
Transfinite Curve{-165} = Ne_y1 Using Progression 1/r_y1;

Line(168) = {18, 68};
Transfinite Curve{-168} = Ne_y1 Using Progression 1/r_y1;

// - Top Vertical lines
Line(209) = {59, 109};
Transfinite Curve{-209} = Ne_y2 Using Progression r_y2;

Line(212) = {62, 112};
Transfinite Curve{-212} = Ne_y2 Using Progression r_y2;

Line(215) = {65, 115};
Transfinite Curve{-215} = Ne_y2 Using Progression r_y2;

Line(218) = {68, 118};
Transfinite Curve{-218} = Ne_y2 Using Progression r_y2;


// - Bottom xz-plane
Line Loop(12) = {13, 26, -20, -23};
Plane Surface(12) = {12};

// - Middle xz-plane
Line Loop(21) = {63, 76, -70, -73};
Plane Surface(21) = {21};

// - Top xz-plane
Line Loop(31) = {113, 126, -120, -123};
Plane Surface(31) = {31};


// - Bottom Vertical Surfaces:
Line Loop(113) = {-13, 159, 63, -162};
Plane Surface(113) = {113};

Line Loop(120) = {-20, 168, 70, -165};
Plane Surface(120) = {120};


Line Loop(123) = {23, 168, -73, -159};
Plane Surface(123) = {123};


Line Loop(126) = {26, 165, -76, -162};
Plane Surface(126) = {126};


// - Top Vertical Surfaces:
Line Loop(114) = {-63, 209, 113, -212};
Plane Surface(114) = {114};

Line Loop(121) = {-70, 218, 120, -215};
Plane Surface(121) = {121};


Line Loop(124) = {73, 218, -123, -209};
Plane Surface(124) = {124};


Line Loop(127) = {76, 215, -126, -212};
Plane Surface(127) = {127};




// - Solids
// - - Bottom:
Surface Loop(12) = {12, 21, 113, 120, 123, 126};
Volume(12) = {12};
// - - Top:
Surface Loop(102) = {21, 31, 114, 121, 124, 127};
Volume(102) = {102};




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
Physical Surface("inlet", 1) = {113,114};
Physical Surface("outlet", 2) = {120,121};
Physical Surface("bottom", 3) = {12};
Physical Surface("top", 4) = {31};
Physical Surface("left", 5) = {123,124};
Physical Surface("right", 6) = {126,127};
Physical Volume("fluid", 8) = {12,102};