############################################
#Author: Andreas Slyngstad
#MEK
#Solving Poission Equation with both Dirichlet
#and Neumann conditions
#############################################

from dolfin import *
import numpy as np
#http://www.math.rutgers.edu/~falk/math575/Boundary-conditions.html



kl = [1, 10, 100]
h = [8, 16, 32, 64]

y = np.zeros(len(h)); x = np.zeros(len(h))
for j in range(len(h)):

    for i in kl:
        mesh = UnitSquareMesh(h[j], h[j])

        #Defining spaces and functions
        V = FunctionSpace(mesh, 'CG', 1)
        u = TrialFunction(V)
        v = TestFunction(V)

        class Dirichlet(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and(near(x[0], 0) or near(x[0], 1) )

        diri = Dirichlet()
        #Setting boundary values
        boundaries = FacetFunction("size_t", mesh)
        boundaries.set_all(0)

        diri.mark(boundaries,1)

        #plot(boundaries); interactive()

        bc0 = DirichletBC(V, 0, boundaries, 1)

        #Defining variational problem
        k = 1; l = 1
        u_e = interpolate(Expression('sin(pi*k*x[0])*cos(pi*l*x[1])', k=k, l=l), V)
        #test = Expression("1")
        a = inner(grad(u), grad(v))*dx
        L = - div( grad(u_e) ) * v * dx   #-nabla(u_e) = f
        #L = u_e*v*dx

        u_ = Function(V)
        solve(a == L, u_, bc0)


        #Norms of the error
        L2 = errornorm(u_e, u_, norm_type='l2', degree_rise = 3)
        H1 = errornorm(u_e, u_, norm_type='H1', degree_rise = 3)

        print "----------------------------------"
        print "For %d points and k, l = %d" % (j, i)
        print "L2 Norm = %.5f -----  H1 Norm = %.5f" % (L2, H1)
        print
        if i == 1:
            d = mesh.coordinates()
            y[j] = np.log(L2 / norm(u_e, 'H1') )
            x[j] = np.log(d[1][0]-d[0][0])


import matplotlib.pyplot as plt
#Error estimates
A = np.zeros((2, 2))
b = np.zeros(2)
A[0][0] = len(h)
A[0][1] = np.sum(x); A[1][0] = A[0][1]
A[1][1] = A[0][1]*A[0][1]
b[0] = np.sum(y); b[1] = np.sum(y*x)

a, b = np.linalg.solve(A, b)
print a

plt.figure(1)
plt.plot(a*x + b)
#plt.plot(y)
plt.show()
