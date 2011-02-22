// mesh width associated with points
lc = 0.06;

// vertices of the cube
Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {0, 1, 0, lc};
Point(4) = {1, 1, 0, lc};
Point(5) = {0, 0, 1, lc};
Point(6) = {1, 0, 1, lc};
Point(7) = {0, 1, 1, lc};
Point(8) = {1, 1, 1, lc};

// lines of the cube
Line(1) = {1,2};
Line(2) = {2,4};
Line(3) = {4,3};
Line(4) = {3,1};
Line(5) = {5,6};
Line(6) = {6,8};
Line(7) = {8,7};
Line(8) = {7,5};
Line(9) = {1,5};
Line(10) = {2,6};
Line(11) = {3,7};
Line(12) = {4,8};

// faces of the cube
Line        Loop(100) = {1,2,3,4};  
Plane    Surface(200) = {100};  
Physical Surface(300) = {200};

Line        Loop(101) = {9,5,-10,-1};  
Plane    Surface(201) = {101};  
Physical Surface(301) = {201};

Line        Loop(102) = {10,6,-12,-2};  
Plane    Surface(202) = {102};  
Physical Surface(302) = {202};

Line        Loop(103) = {12,7,-11,-3};  
Plane    Surface(203) = {103};  
Physical Surface(303) = {203};

Line        Loop(104) = {11,8,-9,-4};  
Plane    Surface(204) = {104};  
Physical Surface(304) = {204};

Line        Loop(105) = {-8,-7,-6,-5};  
Plane    Surface(205) = {105};  
Physical Surface(305) = {205};

// volume of the cube
Surface Loop(400) = {-200,-201,-202,-203,-204,-205};
Volume(1000) = {400}; 
Physical Volume(1001) = {1000};
