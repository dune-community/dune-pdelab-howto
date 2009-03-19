// mesh width associated with points
lc = 1;

// vertices of the pyramid
Point(1) = {-1, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 0, 0, lc};

Point(4) = {2,1,0,lc};
Point(5) = {2.5,3,0,lc};
Point(6) = {-1,1,0,lc};


Circle(1) = {1,3,2};
BSpline(2) = {2,4,5,6,1};

Line Loop(100) = {1,2};  
Plane Surface(200) = {100};  
