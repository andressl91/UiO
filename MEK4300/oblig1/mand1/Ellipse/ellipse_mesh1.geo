cl = 0.2;
Point(2) = {0, 0, 0, cl};
	Point(3) = {-2, 0, 0, cl};
	Point(4) = {2, 0, 0, cl};
	Point(5) = {0, 1, 0, cl};
	Point(6) = {0, -1, 0, cl};
	Ellipse(1) = {6, 2, 4, 4};
	Ellipse(2) = {4, 2, 5, 5};
	Ellipse(3) = {5, 2, 3, 3};
	Ellipse(4) = {3, 2, 6, 6};
	Line Loop(6) = {3, 4, 1, 2};
	Plane Surface(6) = {6};