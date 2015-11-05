from dolfin import *

#Define space and mesh
l = 5.0
h = 1.0
#mesh = RectangleMesh(0.0,0.0,l,h,50,10)
mesh = Mesh("box1.xml")

#Different P1 and P2 due to stability
V = VectorFunctionSpace(mesh, "Lagrange", 2)
Q = FunctionSpace(mesh, "Lagrange", 1)
W = V*Q

#Defining boundaries
class Left(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0],-6.0)

class Right(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[0],6.0)

class Top(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1],2.0)

class Bottom(SubDomain):
	def inside(self, x, on_boundary):
		return near(x[1],-2.0)

class Circle(SubDomain):
    def inside(self, x, on_boundary):
        return on_boundary and -1.8 < x[1] < 1.8 and -5.5 < x[0] < 5.5 
        #return (between(x[1], (-1.8, 1.8)) and between(x[0], (-5.5, 5.5)) )
        #This takes also all points, not just boundary...


left = Left(); right = Right()
top = Top(); bottom = Bottom()
circle = Circle()

#Interior domain Circle
#size_t  creates mesh function of given dimention
domains = FacetFunction("size_t", mesh)
domains.set_all(0)
circle.mark(domains, 1)

#Boundary domain
boundaries = FacetFunction("size_t", mesh)
boundaries.set_all(0)
left.mark(boundaries, 1)
right.mark(boundaries, 2)
top.mark(boundaries, 3)
bottom.mark(boundaries, 4)

bctop = DirichletBC(V, Constant((0,0)), boundaries, 3)
bcbottom = DirichletBC(V, Constant((0,0)), boundaries, 4)
bccirle = DirichletBC(V, Constant((0,0)), domains, 1)
bcleft = DirichletBC(V, Constant((1.0,0.0)), domains, 1)

bcs = [bctop, bcbottom, bccirle, bcleft]

#Define Variational problem
mu = 100
u, p = TrialFunctions(W)
v, q = TestFunctions(W)

a = mu*inner(grad(v),grad(u))*dx - inner(div(v), p)*dx - inner(q, div(u))*dx
#a = mu*inner(grad(v),grad(u))*dx - inner(q, div(u))*dx
up = Function(W)
solve(lhs(a) == rhs(a), up, bcs)

u_, p_ = up.split()
u, v = u_.split()

plot(u_)
interactive()
file = File("stokes.pvd")
file << u


#plot(u_)
#interactive()

#The solve