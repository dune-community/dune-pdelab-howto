// Gmsh project created on Tue Apr  5 23:43:56 2011
Point(1) = {0.4, 0.4, 0.4};
Point(2) = {0.6, 0.4, 0.4};
Point(3) = {0.4, 0.6, 0.4};
Point(4) = {0.4, 0.4, 0.6};
Point(5) = {0.6, 0.6, 0.4};
Point(6) = {0.6, 0.4, 0.6};
Point(7) = {0.4, 0.6, 0.6};
Point(8) = {0.6, 0.6, 0.6};
Line(1) = {1, 2};
Line(2) = {1, 4};
Line(3) = {4, 6};
Line(4) = {6, 2};
Line(5) = {2, 5};
Line(6) = {5, 8};
Line(7) = {8, 6};
Line(8) = {5, 3};
Line(9) = {3, 7};
Line(10) = {7, 8};
Line(11) = {7, 4};
Line(12) = {3, 1};
Line Loop(13) = {2, -11, -9, 12};
Plane Surface(14) = {13};
Line Loop(15) = {8, 9, 10, -6};
Plane Surface(16) = {15};
Line Loop(17) = {5, 6, 7, 4};
Plane Surface(18) = {17};
Line Loop(19) = {3, -7, -10, 11};
Plane Surface(20) = {19};
Line Loop(21) = {4, -1, 2, 3};
Plane Surface(22) = {21};
Line Loop(23) = {12, 1, 5, 8};
Plane Surface(24) = {23};
Point(9) = {0.5, 0.5, 0.5};
Point(10) = {0.0, 0.5, 0.5};
Point(11) = {0.5, 0.0, 0.5};
Point(12) = {0.5, 0.5, 1.0};
Point(13) = {0.5, 0.5, 0.0};
Circle(25) = {10, 9, 11};
Circle(26) = {10, 9, 12};
Circle(27) = {11, 9, 12};
Line Loop(28) = {25, 27, -26};
Ruled Surface(29) = {28};
Rotate {{0, 0, 1}, {0.5, 0.5, 0.5}, Pi/2} {
  Duplicata { Surface{29}; }
}
Rotate {{0, 0, 1}, {0.5, 0.5, 0.5}, Pi/2} {
  Duplicata { Surface{30}; }
}
Rotate {{0, 0, 1}, {0.5, 0.5, 0.5}, Pi/2} {
  Duplicata { Surface{33}; }
}
Circle(38) = {13, 9, 10};
Circle(39) = {13, 9, 11};
Line Loop(40) = {39, -25, -38};
Ruled Surface(41) = {40};
Rotate {{0, 0, 1}, {0.5, 0.5, 0.5}, Pi/2} {
  Duplicata { Surface{41}; }
}
Rotate {{0, 0, 1}, {0.5, 0.5, 0.5}, Pi/2} {
  Duplicata { Surface{42}; }
}
Rotate {{0, 0, 1}, {0.5, 0.5, 0.5}, Pi/2} {
  Duplicata { Surface{44}; }
}
Surface Loop(47) = {46, 41, 42, 44, 33, 36, 29, 30};
Surface Loop(48) = {16, 24, 14, 22, 18, 20};
Volume(49) = {47, 48};
Volume(50) = {48};
