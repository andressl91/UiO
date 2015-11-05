from dolfin import *
h = 1
N = 10

mesh = IntervalMesh(N,0,h)

#Define Functionspace
V = FunctionSpace(mesh, 'Lagrange',1)

#Define boundary condidtions
u_top = Constant(1.)
u_bottom = Constant(0.)
def bottom(x, on_boundary):
	tol = 1E-15
	return on_boundary and \
			abs(x[0] < tol)
def top(x, on_boundary):
	tol = 1E-15
	return on_boundary and \
			abs(tol < (1-x[0]))

bc1 = DirichletBC(V, u_top, top)
bc2 = DirichletBC(V, u_bottom, bottom)

#Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
a = inner(nabla_grad(u), nabla_grad(v))*dx
L = Constant(5)*v*dx

#Compute solution
u_ = Function(V)
solve(a == L, u_, [bc1,bc2])

plot(u_)
interactive()