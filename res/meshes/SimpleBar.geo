// Gmsh project created on Wed Apr 26 16:05:48 2023
unit = 1e0;

// Points
Point(1) = {-1.0*unit, -1.0*unit, 0.0, 0.2*unit};
Point(2) = { 1.0*unit, -1.0*unit, 0.0, 0.2*unit};
Point(3) = { 1.0*unit, 1.0*unit, 0.0, 0.2*unit};
Point(4) = {-1.0*unit, 1.0*unit, 0.0, 0.2*unit};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

// Bar
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Boundaries
Physical Curve("Right", 5) = {2};
Physical Curve("Left", 6) = {4};
Physical Curve("Top", 7) = {3};
Physical Curve("Bottom", 8) = {1};

Physical Surface("Plane", 9) = {1};

// Mesh Definition
Transfinite Surface {1} = {4, 3, 2, 1};
Transfinite Curve {4, 2} = 50 Using Progression 1;
Transfinite Curve {3, 1} = 50 Using Progression 1;

Recombine Surface {1};
//+
Physical Curve("TEST", 10) = {3, 2};
