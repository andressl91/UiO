Point(1) = {0,0,0,lc};
Point(2) = {1,0,0,lc};
Point(3) = {0,1,0,lc};
Point(4) = {-1,0,0,lc};

Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Line(3) = {4,2};

Line Loop(4) = {1,2,3}; 
Plane Surface(5) = {4};
