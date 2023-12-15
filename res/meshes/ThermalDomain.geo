// Gmsh project created on Sat May 20 12:32:18 2023
// Dimensions and mesh size
units = 12.5e-3;
he = 1.0 * units;


// Points
Point(1) = { 0 * units, 0 * units, 0 * units, he };
Point(2) = { 15 * units, 0 * units, 0 * units, he };
Point(3) = { 17 * units, 0 * units, 0 * units, he };
Point(4) = { 32 * units, 0 * units, 0 * units, he };
Point(5) = { 32* units, 1 * units, 0 * units, he };
Point(6) = { 32 * units, 5 * units, 0 * units, he };
Point(7) = { 32 * units, 16 * units, 0 * units, he };
Point(8) = { 0 * units, 16 * units, 0 * units, he };

// Edges
Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 5 };
Line(5) = { 5, 6 };
Line(6) = { 6, 7 };
Line(7) = { 7, 8 };
Line(8) = { 8, 1 };

// Domain Definition
Curve Loop(1) = { 1, 2, 3, 4, 5, 6, 7, 8 };
Plane Surface(1) = { 1 };

// Boundaries
Physical Point("SouthWestNode", 1) = {1};
Physical Point("SouthEastNode", 2) = {4};

Physical Curve("Heat", 3) = {2};
Physical Curve("Left", 4) = {8};
Physical Curve("Right", 5) = {6, 5, 4};
Physical Curve("Top", 6) = {7};
Physical Curve("Bottom", 7) = {1, 2, 3};
Physical Curve("TractionBottom", 8) = {4};
Physical Curve("TractionTop", 9) = {6};

Physical Surface("Domain", 10) = { 1 };

// Mesh
//Recombine Surface {1};