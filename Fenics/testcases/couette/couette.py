author = "Andreas Slyngstad"
date = "11/01/2016"

from dolfin import *

mesh = IntervalMesh(50, 0, 2)
V = FunctionSpace(mesh, 'CG', 1)

u = TrialFunction(V)
v = TestFunction(V)

def bottom(x, on_boundary):
	return on_boundary and abs(x[0] < 1E-8)
def top(x, on_boundary):
	return on_boundary and abs((1-x[0]) < 1E-8)

bc1 = DirichletBC(V, Constant(1), top)
bc2 = DirichletBC(V, Constant(0), bottom)
bcs = [bc1, bc2]

#Variational Form
a = inner(grad(u), grad(v))*dx
L = Constant(0)*v*dx

u_ = Function(V)
solve(a == L, u_, bcs)

plot(u_)
interactive()