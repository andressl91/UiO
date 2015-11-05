from dolfin import *

#Define space and mesh
l = 5.0
h = 1.0
mesh = RectangleMesh(0.0,0.0,l,h,50,10)
#plot(mesh)
#interactive()

V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
W = V*Q

#Define Boundary Conditions
def top(x, on_boundary):
	tol = 1E-15
	return on_boundary or \
			abs(tol < 1-x[1])
			
def bottom(x, on_boundary):
	tol = 1E-15
	return on_boundary or \
			abs(tol < x[1])

def left(x, on_boundary):
	tol = 1E-15
	return on_boundary or \
			abs(tol < x[0])

def right(x, on_boundary):
	tol = 1E-15
	return on_boundary or \
			abs(tol < 1-x[0])

btop = DirichletBC(V, 1.0, top)
bbottom = DirichletBC(V, 0.0, bottom)
bright

#Define Variational problem
#up = Function(W)
#(u,p) = split(up)

u = Function(W) 
p = Function(W)
(u, p) = TrialFunctions(W)
(v, q) = TestFunctions(W)
