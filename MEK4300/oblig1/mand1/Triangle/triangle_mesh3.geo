cl = 0.0666666666667;
Point(1) = {0, 0, 0, cl};
				Point(2) = {-0.5, Sqrt(0.75), 0, cl};
				Point(3) = {0.5, Sqrt(0.75), 0, cl};
				Line(1) = {1, 3};
				Line(2) = {3, 2};
				Line(3) = {2, 1};
				Line Loop(4) = {3, 1, 2};
				Plane Surface(5) = {4};