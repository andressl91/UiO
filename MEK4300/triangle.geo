cl = 0.1
Point(1) = {0,0,0,cl};
Point(2) = {0,0.5,0.866025,cl};
Point(3) = {0,-0.5,0.866025,cl};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,1};

Line Loop(4) = {1,2,3};
Plane Surface(5) = {4};
