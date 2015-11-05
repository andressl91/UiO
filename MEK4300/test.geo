lc = 1.0/128.0;

Point(1) = {0.0, 0.0, 0, lc};
Point(2) = {0.1, 0.0, 0, lc};
Point(3) = {0.5, -0.1,  0, lc} ; //Ellipse center point
Point(4) = {0.5, 0.1, 0, lc} ;  // Ellipse top point
Point(5) = {0.9, 0.0, 0, lc};
Point(6) = {1.0,  0.0, 0, lc} ;
Point(7) = {1.0, 0.3, 0, lc};
Point(8) = {0.9,  0.3, 0, lc} ;
Point(9) = {0.5, 0.2,  0, lc} ;  //Ellipse top point
Point(10) = {0.5, 0.4, 0, lc} ;  //Ellipse center point
Point(11) = {0.1, 0.3, 0, lc};
Point(12) = {0.0, 0.3, 0, lc} ;

Line(1) = {1,2};
Ellipse(2) = {4,3,3,2} ;
Ellipse(3)={4,3,3,5} ;
Line(4) = {5,6} ;
Line(5) = {6,7};
Line(6) = {7,8};
Ellipse(7) = {9,10,10,8} ;
Ellipse(8) = {9,10,10,11} ;
Line(9) = {11,12};
Line(10) = {12,1} ;

Line Loop(11) = {1,-2,3,4,5,6,-7,8,9,10} ;

Plane Surface(12) = {11} ;
