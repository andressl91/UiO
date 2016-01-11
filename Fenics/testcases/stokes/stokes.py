from dolfin import *
import numpy as np
#Stokes challenge http://math.tifrbng.res.in/~fem2012/FeNICS_files/tifr_lecture_05.pdf

mesh = Mesh('dolfin_fine.xml')

V = VectorFunctionSpace(mesh, 'Lagrange',2)
Q = FunctionSpace(mesh,'Lagrange',1)
W = V*Q

class Left(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[0],0)

class Right(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[0],1)

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[1],1)

class Bottom(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and near(x[1],0)

class Dolfin(SubDomain):
	def inside(self, x, on_boundary):
		return on_boundary and not near(x[0],0) \
				and not near(x[0],1) \
				and not near(x[1],0) \
				and not near(x[1],1)

left = Left(); right = Right()
top = Top(); bottom = Bottom()
dolf = Dolfin()

#Boundary domain
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left.mark(boundaries, 1)
right.mark(boundaries, 2)
top.mark(boundaries, 3)
bottom.mark(boundaries, 4)
dolf.mark(boundaries,5)
#plot(boundaries)
#interactive()
from math import pi
u_in = interpolate(Expression(('-sin(pi*x[1])','0'), pi=pi), V)
p_in = interpolate(Expression('0'), Q)


bcleft = DirichletBC(Q, p_in, boundaries, 1)

bcright = DirichletBC(V, u_in, boundaries, 2)
bctop = DirichletBC(V, Constant((0,0)), boundaries, 3)
bcbottom = DirichletBC(V, Constant((0,0)), boundaries, 4)
bcdolfin = DirichletBC(V, Constant((0,0)), boundaries, 5)

bcs = [bctop, bcbottom, bcright, bcleft, bcdolfin]

nu = 1000
u, p = TrialFunctions(W)
v, q = TestFunctions(W)

#Variational Form
a = inner(grad(v),nu*grad(u))*dx - p*div(v)*dx - q*div(u)*dx

u_p = Function(W)
solve(lhs(a) == rhs(a), u_p, bcs)
u_, p_ = split(u_p)

plot(u_)
interactive()
