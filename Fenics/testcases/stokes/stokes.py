from dolfin import *

#SJEKK DENNE http://www.jstor.org/stable/2100799?seq=1#page_scan_tab_contents
 #Define space and mesh
l = 5.0
h = 1.0
#mesh = RectangleMesh(0.0,0.0,l,h,50,10)
mesh = Mesh("coin.xml")

#Different P1 and P2 due to stability
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
W = V*Q

#Defining boundaries
class Left(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0],0.0)

class Right(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0],1.0)

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1],1.0)

class Bottom(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1],0.0)

class Dolf(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and not near(x[0],0.0) \
        					and not near(x[0],1.0) \
        					and not near(x[1],0.0) \
        					and not near(x[1],1.0)


left = Left(); right = Right()
top = Top(); bottom = Bottom()
dolfin = Dolf()

#Interior domain Circle

#Boundary domain
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left.mark(boundaries, 1)
right.mark(boundaries, 2)
top.mark(boundaries, 3)
bottom.mark(boundaries, 3)
dolfin.mark(boundaries, 3)

bcnoslip = DirichletBC(W.sub(0), Constant((0,0)), boundaries, 3)
bcright = DirichletBC(W.sub(0), Constant((-5.0,0.0)), boundaries, 2)

#plot(boundaries); interactive()

bcs = [bcnoslip,  bcright]

#Define Variational problem
mu = 0.1
u, p = TrialFunctions(W)
v, q = TestFunctions(W)

a = mu*inner(grad(v),grad(u))*dx - inner(div(v), p)*dx + inner(q, div(u))*dx
#a = mu*inner(grad(v),grad(u))*dx - inner(q, div(u))*dx
up = Function(W)
solve(lhs(a) == rhs(a), up, bcs)

u_, p_ = up.split()
u, v = u_.split()

file = File("stokes.pvd")
file << u_

plot(u_)
interactive()
