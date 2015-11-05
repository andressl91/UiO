from dolfin import *
import os
import sys
from math import *

#Setting Constants
mu = 2.
dpdx = -0.1
a = 2.
b = 1.

#Defining exact solution

y = Expression('x[1]')
x = Expression('x[0]')
u_e = 1.0/(2.*mu)*(-dpdx)*((a*a*b*b)/(a*a + b*b))*(1-x*x/(a*a) - y*y/(b*b))

#Defining lists 
error = [] #Error
h = []
cls = [] 

#Function to make ellipses with a certain mesh size
def makeel(filen, cl):
	file = open(filen, "w")

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

#IMPORTANT n is NUMBER OF MESHES MADE, FIXED TO THE FUNCTIONS AND LOOPS

def divide(cl, n):
	#Reducing the mesh size cl=0.5 with a half, n times
	#Remember that the first calculation is NOT cl, but cl divided in two

	for i in range(1,n+1): #n+1 to get the right amount of files to make 
		calc = cl/(2.0*i)
		cls.append(calc)		
		inp = str(i)
		filename = 'ellipse_mesh%s' % inp
		makeel('%s.geo' % filename, calc)
		os.system('gmsh -2 %s.geo' % filename)
		os.system('dolfin-convert %s.msh %s.xml' % (filename, filename))


def main(deg, mesh):
	#Define variational problems
	mesh = mesh
	V = FunctionSpace(mesh, 'CG', deg)
	u = TrialFunction(V)
	v = TestFunction(V)
	F = inner(grad(u),grad(v))*dx + 1.0/mu*dpdx*v*dx
	bc = DirichletBC(V, Constant(0), DomainBoundary()) 

	#Solving the problem
	u_ = Function(V)
	solve(lhs(F) == rhs(F), u_, bc)

	#Calculating Error
	u_er = project(u_e, V)
	err = errornorm(u_er, u_, degree_rise=0)
	error.append(err)
	interactive()


def calc(n, highd):
	for k in range(1,n+1):
		mesh = Mesh('ellipse_mesh%d.xml' % k) 	#Currrent mesh size to calculate
		hi = mesh.hmin()
		h.append(hi)
		for j in range(1, highd+1): 		  	#Solving for lagrange from degree 1 up to highp
			main(j, mesh)					    #Appending the errors in the lists


def compare(highp, n):
	for j in range(n-1):		
		for i in range(highp):	
				r = log( error[i+highp*(j+1)] / error[i+highp*(j)] ) / log( h[j+1] / h[j] )
				
				print "h=%2.2E Error=%2.2E r=%.2f" % (h[j+1], error[i+highp*(j+1)], r)

n=4				    #Number of times to divide cl in half
cl = 0.4		    #Mesh point size
deg = 3				#Highest Lagrange degree to be calculated

divide(cl, n)       #Read as divide cl, n times
calc(n, deg)	    #Read as calcualte n meshes, with lagrange degree up to deg 	
compare(deg, n) 	#Read as compare the results with lagrange degree from 1 to deg, with the corresponding meshes n




