from dolfin import *

mesh = Mesh("box1.xml")
plot(mesh)
interactive()
#Define function spaces
V = VectorFunctionSpace(mesh, "CG", 2)
Q = FunctionSpace(mesh, "CG", 1)

#Define Trial and Test Function
u = TrialFunction(V)
p = TrialFunction(Q)
v = TestFunction(V)
q = TestFunction(Q)

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


left = Left()
right = Right()
top = Top()
bottom = Bottom()
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

#Set paramters
dt = 0.01
T = 10
nu = 0.01

# Define time-dependent pressure boundary condition
p_in = Expression("1.0", t=0.0)


bctop = DirichletBC(V, Constant((0.0,0.0)), boundaries, 3)
bcbottom = DirichletBC(V, Constant((0.0,0.0)), boundaries, 4)
bccirle = DirichletBC(V, Constant((0.0,0.0)), domains, 1)
bcleft = DirichletBC(V, Constant((1.0,0.0)), domains, 1)

inflow  = DirichletBC(Q, p_in, boundaries, 1)
outflow = DirichletBC(Q, p_in, boundaries, 2)

bcu = [bctop, bcbottom, bccirle, bcleft]
bcp = [inflow, outflow]

# Create functions
u0 = Function(V)
u1 = Function(V)
p1 = Function(Q)

# Define coefficients
k = Constant(dt)
f = Constant((0, 0))

# Tentative velocity step
F1 = (1/k)*inner(u - u0, v)*dx + inner(grad(u0)*u0, v)*dx + \
     nu*inner(grad(u), grad(v))*dx - inner(f, v)*dx
a1 = lhs(F1)
L1 = rhs(F1)

# Pressure update
a2 = inner(grad(p), grad(q))*dx
L2 = -(1/k)*div(u1)*q*dx

# Velocity update
a3 = inner(u, v)*dx
L3 = inner(u1, v)*dx - k*inner(grad(p1), v)*dx

# Assemble matrices
A1 = assemble(a1)
A2 = assemble(a2)
A3 = assemble(a3)

# Use amg preconditioner if available
prec = "amg" if has_krylov_solver_preconditioner("amg") else "default"

# Create files for storing solution
ufile = File("results/velocity.pvd")
pfile = File("results/pressure.pvd")

# Time-stepping
t = dt
while t < T + DOLFIN_EPS:

    # Update pressure boundary condition
    p_in.t = t

    # Compute tentative velocity step
    begin("Computing tentative velocity")
    b1 = assemble(L1)
    [bc.apply(A1, b1) for bc in bcu]
    solve(A1, u1.vector(), b1, "gmres", "default")
    end()

    # Pressure correction
    begin("Computing pressure correction")
    b2 = assemble(L2)
    [bc.apply(A2, b2) for bc in bcp]
    solve(A2, p1.vector(), b2, "gmres", prec)
    end()

    # Velocity correction
    begin("Computing velocity correction")
    b3 = assemble(L3)
    [bc.apply(A3, b3) for bc in bcu]
    solve(A3, u1.vector(), b3, "gmres", "default")
    end()

    # Plot solution
    #plot(p1, title="Pressure", rescale=True)
    #plot(u1, title="Velocity", rescale=True)

    # Save to file
    ufile << u1
    pfile << p1

    # Move to next time step
    u0.assign(u1)
    t += dt
    print "t =", t

# Hold plot
#interactive()