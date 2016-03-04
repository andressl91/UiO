############################################
#Author: Andreas Slyngstad
#MEK
#Solving Poission Equation with both Dirichlet
#and Neumann conditions
#############################################

from dolfin import *
#http://www.math.rutgers.edu/~falk/math575/Boundary-conditions.html

mesh = UnitSquareMesh(10, 10)

#Defining spaces and functions
V = FunctionSpace(mesh, 'CG', 2)
u = TrialFunction(V)
v = TestFunction(V)

def Dirichlet(x, on_boundary):
    return on_boundary and( near(x[0],0) or near(x[0], 1) )

bc0 = DirichletBC(V, 0, Dirichlet)
#Defining variational problem
k= 1; l=1
u_e = Expression('sin(pi*k*x[0])*cos(pi*l*x[1])', k=k, l=l)
a = inner(grad(u), grad(v))*dx
L = -div(grad(u_e))*v*dx #-nabla(u_e) = f

u_ = Function(V)
solve(a == L, u_, bc0)

#H^P norm of u
#H_p = sqrt(assemble( ) )

#Norms of the error
err = u_e - u
