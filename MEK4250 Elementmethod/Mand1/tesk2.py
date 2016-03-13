############################################
#Author: Andreas Slyngstad
#MEK
#Solving Poission Equation with both Dirichlet
#and Neumann conditions
#############################################

from dolfin import *
import numpy as np
#http://www.math.rutgers.edu/~falk/math575/Boundary-conditions.html


class Poission():
    def __init__(self, h):
        self.y = np.zeros(len(h))
        self.x = np.zeros(len(h))
        self.count = 0

    def set_mesh(self,i):
        self.h = i
        self.mesh = UnitSquareMesh(i, i)

    def calc(self, k, l, output=True):
        mesh = self.mesh

        #Defining spaces and functions
        V = FunctionSpace(mesh, 'CG', 1)
        u = TrialFunction(V)
        v = TestFunction(V)

        class Bottom(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], 0)
        class Top(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], 1)

        top = Top()
        bottom = Bottom()
        #Setting boundary values
        boundaries = FacetFunction("size_t", mesh)
        boundaries.set_all(0)
        bottom.mark(boundaries,1)
        top.mark(boundaries, 2)
        bc0 = DirichletBC(V, 0, bottom)
        bc1 = DirichletBC(V, 1, top)

        #Defining and solving variational problem
        u_e = interpolate(Expression('sin(pi*k*x[0])*cos(pi*l*x[1])', k=k , l=l), V)
        f = Expression("((pi*pi*k*k)+(pi*pi*l*l))*sin(pi*k*x[0])*cos(pi*l*x[1])",k=k,l=l)
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx
        #L = -div( grad(u_e))*v*dx   #-nabla(u_e) = f

        u_ = Function(V)
        solve(a == L, u_, bc0, bc1)

        #Norms of the error
        L2 = errornorm(u_e, u_, norm_type='L2', degree_rise = 3)
        H1 = errornorm(u_e, u_, norm_type='H1', degree_rise = 3)
        if output == True:
            print "----------------------------------"
            print "For %d points and k, l = %d" % (self.h, k)
            print "L2 Norm = %.5f -----  H1 Norm = %.5f" % (L2, H1)
            print
        if k == 1:
            d = mesh.coordinates()
            self.x[self.count] = np.log(1./self.h)
            self.y[self.count] = np.log( L2 / norm(u_e, "H1") )
            self.count += 1

    def l_square(self, h, fig):
        A = np.zeros((2, 2))
        b = np.zeros(2)

        A[0][0] = len(h)
        A[0][1] = np.sum(self.x); A[1][0] = A[0][1]
        A[1][1] = np.sum(self.x*self.x)
        b[0] = np.sum(self.y)
        b[1] = np.sum(self.y*self.x)

        a, b = np.linalg.solve(A, b)

        if fig == True:
            import matplotlib.pyplot as plt
            plt.figure(1)
            plt.plot(self.x, b*self.x + a, label='Linear approximation')
            plt.plot(self.x, self.y, 'o', label='Points to be approximated')
            plt.legend()
            plt.show()

kl = [1, 10, 100]
h = [2**i for i in range(4)]

prob = Poission(h)
for i in h:
    for k in kl:
        prob.set_mesh(i)
        prob.calc(k,k, output = True)
#prob.l_square(h, fig = True)
