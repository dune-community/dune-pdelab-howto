// mesh width associated with points
lc = 0.05;

//          5             4
//     6 ---5-------- 4
//     |           \              |
//     |             7 \         | 3
//     |                    \     |
// 6  |                           3
//     |                           |
//     |                           | 2
//     |                           |
//     1 ------------ 2
//                     1

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 0.4, 0, lc};
Point(4) = {1, 1, 0, lc};
Point(5) = {0.1, 1.0, 0, lc};
Point(6) = {0.0, 1.0, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,5};
Line(5) = {5,6};
Line(6) = {6,1};
Line(7) = {3,5};

Line Loop(100) = {1,2,7,5,6};  
Line Loop(101) = {3,4,-7};  

Plane Surface(200) = {100};  
Plane Surface(201) = {101};  

