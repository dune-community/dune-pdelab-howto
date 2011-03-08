// mesh width associated with points
lc = 0.1;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {-1, 1, 0, lc};
Point(5) = {-1, -1, 0, lc};
Point(6) = {0, -1, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};

Line Loop(100) = {1,2,3,4,5,6};  
Plane Surface(200) = {100};  
