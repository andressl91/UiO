from dolfin import *
import numpy


T = 10.0 # tenison
A = 1.0  # pressure amplitude
R = 0.3  # radius of domain
theta = 0.2
x0 = 0.6*R*cos(theta)
yo = 0.6*R*sin(theta)
sigma = 0.025
n = 40

# Making the mesh
mesh = UnitCircleMesh(n) #n denotes elements along radial dir
V = FunctionSpace(mesh, 'Lagrange', 1)

#Define the BCD w = 0
def boundary(x, on_boundary):
	return on_boundary

bc = DirichletBC(V, Constant(0,0), boundary)
#Sets DiricletBCD, sets the boundary values = 0

w = TrialFunction(V)
v = TestFunction(V)
a = inner(nabla_grad(w), nabla_grad(v))*dx

f = Expression('4*exp(-0.5*(pow((R*x[0] - x0)/sigma, 2)) '
               '     - 0.5*(pow((R*x[1] - y0)/sigma, 2)))',
               R=R, x0=x0, y0=y0, sigma=sigma)
L = f*v*dx
#pow denotes the exponent used on the base


# Compute solution
w = Function(V)
problem = LinearVariationalProblem(a, L, w, bc)
solver  = LinearVariationalSolver(problem)
solver.parameters['linear_solver'] = 'cg'
solver.parameters['preconditioner'] = 'ilu'
solver.solve()


# Plot scaled solution, mesh and pressure
plot(mesh, title='Mesh over scaled domain')
plot(w, title='Scaled deflection')
f = interpolate(f, V)
plot(f, title='Scaled pressure')

interactive()

