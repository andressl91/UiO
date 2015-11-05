from dolfin import *

mu = 100
mesh = UnitSquareMesh(20,20)

#Different degree, due to unstabiity if the degree is the same
V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)

VQ = V*Q
u, p = TrialFunctions(VQ)
v, q = TestFunctions(VQ)

F = mu*inner(grad(v),grad(u))*dx - inner(div(v), p)*dx - inner(q, div(u))*dx
#Negative velocity divergence due to sym coefficient matrix, better iteration
def top(x, on_boundary):
	return x[1] > (1-DOLFIN_EPS)

bc0 = DirichletBC(VQ.sub(0), Constant((0.0,0)), DomainBoundary())
bc1 = DirichletBC(VQ.sub(0), Constant((1.0,0)), top)

up_ = Function(VQ)
solve(lhs(F) == rhs(F), up_, [bc0, bc1])
u_, p_ = up_.split()

u, v = u_.split()

w = v.dx(0) - u.dx(1) #Equation 136
psi = TrialFunction(Q)
psi_v = TestFunction(Q) 

lhs = -inner(grad(psi_v), grad(psi))*dx #Equation 162
rhs = -dot(psi_v, w)*dx

bc = DirichletBC(Q, Constant(0), DomainBoundary())
psi_s = Function(Q)

solve(lhs == rhs, psi_s, bc)

vortex = psi_s.vector().array().argmin() #Lowest point
print(mesh.coordinates()[vortex])

viz = plot(psi_s)
viz.write_png('stokes')

interactive()

viz = plot(u_)
viz.write_png('stokes1')
#interactive()
