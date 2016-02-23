from dolfin import *
import numpy as np
mu = 100.0

mesh = Mesh('figure.xml')

#Different degree, due to unstabiity if the degree is the same
V = VectorFunctionSpace(mesh, 'CG', 2)
Q = FunctionSpace(mesh, 'CG', 1)

VQ = V*Q
u, p = TrialFunctions(VQ)
v, q = TestFunctions(VQ)

Fu = mu*inner(grad(v),grad(u))*dx - inner(div(v), p)*dx - inner(q, div(u))*dx

#Negative velocity divergence due to sym coefficient matrix, better iteration
def top(x, on_boundary):
	return x[1] > (0.5-DOLFIN_EPS) and on_boundary

def bottom(x, on_boundary):
	return (x[1] < (0.1+DOLFIN_EPS) or abs(x[0]-0.5) < DOLFIN_EPS or x[1] < DOLFIN_EPS) and on_boundary

bc0 = DirichletBC(VQ.sub(0), Constant((0,0)), bottom)
bc1 = DirichletBC(VQ.sub(0), Constant((1,0)), top)

up_ = Function(VQ)
solve(lhs(Fu) == rhs(Fu), up_, [bc0, bc1])
u_, p_ = up_.split()

#plot(u_)
#interactive()

u, v = u_.split()

w = v.dx(0) - u.dx(1) #Equation 136

psi = TrialFunction(Q)
psi_v = TestFunction(Q)

#Fetching the normal of the mesh
n = FacetNormal(mesh)
grad_psi = as_vector((-v,u))

#Equation 162 WITH weak boundary
F = -inner(grad(psi_v), grad(psi))*dx + inner(psi_v*grad_psi, n)*ds + psi_v*w*dx

psi_s = Function(Q)

solve(lhs(F) == rhs(F), psi_s)

#viz = plot(psi_s)
#viz.write_png('stokesa')
#interactive()


bc = DirichletBC(Q, Constant(100), DomainBoundary())
bc.apply(psi_s.vector())
val = abs(psi_s.compute_vertex_values())

loc = np.where(val == val.min())
#print mesh.coordinates()[loc]


#Task B
#f = File("psi_s.pvd")
#f << psi_s

#Task C
def left(x, on_boundary):
	return x[0] < 1e-12 and on_boundary
def right(x, on_boundary):
	return x[0] > 1-DOLFIN_EPS and on_boundary

Left = AutoSubDomain(left)
mf = FacetFunction("size_t", mesh)
mf.set_all(0)
Left.mark(mf, 1)

Right = AutoSubDomain(right)
mf.set_all(0)
Right.mark(mf, 1)

Leftintegral = assemble(dot(u_, -n)*ds, exterior_facet_domains=mf)
Rightintegral = assemble(dot(u_, n)*ds, exterior_facet_domains=mf)
#print "The flux on the left boundary is %e" % Leftintegral
#print "The flux on the right boundary is %e" % Rightintegral

#Task D
def reverse():
	bc1 = DirichletBC(VQ.sub(0), Constant((-1,0)), top)

	up_ = Function(VQ)
	solve(lhs(Fu) == rhs(Fu), up_, [bc0, bc1])
	u_, p_ = up_.split()

	u, v = u_.split()

	w = v.dx(0) - u.dx(1) #Equation 136

	psi = TrialFunction(Q)
	psi_v = TestFunction(Q)

	#Fetching the normal of the mesh
	n = FacetNormal(mesh)
	grad_psi = as_vector((-v,u))

	#Equation 162 WITH weak boundary
	F = -inner(grad(psi_v), grad(psi))*dx + inner(psi_v*grad_psi, n)*ds + psi_v*w*dx

	psi_s = Function(Q)

	solve(lhs(F) == rhs(F), psi_s)
	#plot(psi_s)
	#interactive()

	bc = DirichletBC(Q, Constant(100), DomainBoundary())
	bc.apply(psi_s.vector())
	val = abs(psi_s.compute_vertex_values())

	loc = np.where(val == val.min())
	print mesh.coordinates()[loc]

reverse()

#Task E

tau = -p_*Identity(2) + mu*(grad(u_) + grad(u_).T)

Bottom = AutoSubDomain(bottom)
mf.set_all(0)
Bottom.mark(mf, 1)

shear = assemble(dot(dot(tau, n), n)*ds, exterior_facet_domains=mf)
print shear
