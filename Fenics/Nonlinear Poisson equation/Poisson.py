# Solving nonlinear poisson equation
#-del * (q(u)del(u)) = f(u)
from dolfin import * 

#Define mesh and functionspace
mesh = UnitIntervalMesh(10) 
V = FunctionSpace(mesh, "CG", 1)

#Define functions in the equation
def q(u):
	return 1.0+u

#Define boundary conditions
def left(x, on_boundary):
	return x[0] < DOLFIN_EPS
def right(x, on_boundary):
	return x[0] > 1-DOLFIN_EPS

bc0 = DirichletBC(V, Constant(0), left) 
bc1 = DirichletBC(V, Constant(1), right)
bcs = [bc0, bc1]

#Define Trial and Testfunctions
u = TrialFunction(V)
v = TestFunction(V)



method = "Picard"

if method == "Newton":
	u_ = interpolate(Expression("x[0]"), V)
	F = inner(q(u_)*grad(u_), grad(v))*dx #f is set to zero - inner(f(u_), v)*dx
	#After this we can simply put solve(F == 0, u_, bcs), which will solve auto with newton
	#Manually computing the mothod
	#Find the Jacobian
	J = derivative(F, u_, u)
	du = Function(V)
	error = 1; k = 0

	while k < 50 and error > 1e-10:
		A = assemble(J)
		b = assemble(-F)
		[bc.apply(A, b, u_.vector()) for bc in bcs]
		solve(A, du.vector(), b)
		u_.vector().axpy(1., du.vector())
		#Add multiple of given vector (AXPY operation)
		error = norm(du.vector())
		print "Error = ", error
		k += 1

if method == "Picard":
	F = inner(q(u_)*grad(u), grad(v))*dx 
	u_1 = Function(V)
	error = 1; k = 0
	while k < 50 and error > 1e-10:
		solve(F == Constant(0)*v*dx, u_, bcs)
		error = errornorm(u_, u_1)
		u_1.assign(u_)
		print "Error = %.5f after %d Iterations " % (error, k)
		k += 1