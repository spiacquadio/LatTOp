// Gmsh project created on Mon Jun 26 15:54:09 2023
he = 0.1;

Point(1) = { 0.0, 0.0, 0.0, he };
Point(2) = { 1.0, 0.0, 0.0, he };
Point(3) = { 1.0, 1.0, 0.0, he };
Point(4) = { 0.0, 1.0, 0.0, he };

Line(1) = { 1, 2 };
Line(2) = { 2, 3 };
Line(3) = { 3, 4 };
Line(4) = { 4, 1 };

Curve Loop(1) = { 1, 2, 3, 4 };
Plane Surface(1) = { 1 };

Physical Curve("Bottom", 5) = {1};
Physical Curve("Top", 6) = {3};
Physical Curve("Left", 7) = {4};
Physical Curve("Right", 8) = {2};
Physical Surface("Domain", 9) = {1};

Transfinite Surface {1} = {4, 3, 2, 1};
Transfinite Curve {4, 2} = 51 Using Progression 1;
Transfinite Curve {3, 1} = 51 Using Progression 1;
Recombine Surface {1};