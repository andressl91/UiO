Point(1) = {-5, 1, 0, 1.0};
Point(2) = {-5, -1, 0, 1.0};
Point(3) = {5, -1, 0, 1.0};
Point(4) = {-5, -1, 0, 1.0};
Point(5) = {5, -1, 0, 1.0};
Point(6) = {5, 1, 0, 1.0};
Line(1) = {1, 6};
Line(2) = {6, 3};
Line(3) = {3, 2};
Line(4) = {2, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
