from dolfin import *

#Construct mesh and define functionspace
mesh = UnitSquareMesh(6, 4)
V = FunctionSpace(mesh, "Lagrange", 1)

#Define boundary condtitions
u0 = Expression("1+x[0]*x[0]+2*x[1]*x[1]")

def u0_boundary(x,on_boundary):
	return on_boundary

bc = DirichletBC(V, u0, u0_boundary)

#Define variational problem
u = TrialFunction(V)
v = TestFunction(V)
f = Constant(-6.0)
a = inner(nabla_grad(u),nabla_grad(v))*dx
L = f*v*dx

#Compute solution
u = Function(V)
solve(a == L ,u , bc)

#Dolfin array
u_nodal_values = u.vector()

#Numpy array
u_array = u_nodal_values.array()


#Plot function
#plot(u)
#interactive()

#Save solution for paraview sim
#file = File("poisson.pvd")
#file << u

#Checking error
#Interpolate to get functionvalues on V
ue = interpolate(u0,V)
ue_array = ue.vector().array()
#gradu=project(grad(u), VectorFunctionSpace(mesh,"Lagrange", 1)) Also a way
#Normalization
max_u = u_array.max()
u_array = u_array/max_u
#u.vector()[:] = u_array Alternative
u.vector().set_local(u_array)
#print u.vector().array()