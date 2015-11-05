from dolfin import *
import os
import sys
from math import *

#Setting constants, runs smoother in C++
mu = 0.2
dpdx = -1.0
a = 1.0

#Function to make triangle with a certain mesh size
def makeel(filen, cl):
	file = open(filen, "w")

	file.write("cl = %s;\n" % cl) 
	file.write("Point(1) = {0, 0, 0, cl};\n\
				Point(2) = {-0.5, Sqrt(0.75), 0, cl};\n\
				Point(3) = {0.5, Sqrt(0.75), 0, cl};\n\
				Line(1) = {1, 3};\n\
				Line(2) = {3, 2};\n\
				Line(3) = {2, 1};\n\
				Line Loop(4) = {3, 1, 2};\n\
				Plane Surface(5) = {4};")
	file.close

#Analythical solution

u_e = Expression('-dpdx/(2.0*sqrt(3.0)*a*mu)*(x[1]-0.5*a*sqrt\
(3.0))*(3.0*x[0]*x[0]-x[1]*x[1])', mu=mu, a=a, dpdx=dpdx)

#Defining lists 
error = [] #Error
h = []	#
cls = [] #Size of mesh points

#IMPORTANT n is NUMBER OF MESHES MADE, FIXED TO THE FUNCTIONS AND LOOPS

def divide(cl, n):
	#Reducing the mesh size cl=0.5 with a half, n times
	#Remember that the first calculation is NOT cl, but cl 
	#divided in two

	for i in range(1,n+1):
		calc = cl/(2.0*i)
		cls.append(calc)		
		inp = str(i)
		filename = 'triangle_mesh%s' % inp
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

	u_er = project(u_e, V)
	
	#Calculating Error
	u_er = project(u_e, V)
	err = errornorm(u_er, u_, degree_rise=0)
	error.append(err)
	


def calc(n, highd):
	for k in range(1,n+1):
		mesh = Mesh('triangle_mesh%d.xml' % k) 	#Currrent mesh size to calculate
		hi = mesh.hmin()
		h.append(hi)
		for j in range(1, highd+1): 		  	
			main(j, mesh)					   


def compare(highp, n):
	for j in range(n-1):		
		for i in range(highp):	
				r = log( error[i+highp*(j+1)] / error[i+highp*(j)] ) / log( h[j+1] / h[j] )
				
				print "h=%2.2E Error=%2.2E r=%.2f" % (h[j+1], error[i+highp*(j+1)], r)

n=4		    		#Number of times to divide cl in half
cl = 0.4		    #Mesh point size
deg = 3				#Highest Lagrange degree to be calculated

divide(cl, n)       #Read as divide cl, n times
calc(n, deg)	    #Read as calcualte n meshes, with lagrange degree up to deg 	
compare(deg, n) 	#Read as compare the results with lagrange degree from 1 to deg, with the corresponding meshes n




