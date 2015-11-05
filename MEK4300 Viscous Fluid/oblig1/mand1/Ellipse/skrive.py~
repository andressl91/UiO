import sys

cl = sys.argv[1]
file = open("ellipsetest_%s.geo" % cl, "w")

file.write("cl = %s;\n" % cl) 
file.write("Point(2) = {0, 0, 0, cl};\n\
Point(3) = {-2, 0, 0, cl};\n\
Point(4) = {2, 0, 0, cl};\n\
Point(5) = {0, 1, 0, cl};\n\
Point(6) = {0, -1, 0, cl};\n\
Ellipse(1) = {6, 2, 4, 4};\n\
Ellipse(2) = {4, 2, 5, 5};\n\
Ellipse(3) = {5, 2, 3, 3};\n\
Ellipse(4) = {3, 2, 6, 6};\n\
Line Loop(6) = {3, 4, 1, 2};\n\
Plane Surface(6) = {6};")
file.close
