############################################
#Author: Andreas Slyngstad
#MEK
#Solving Poission Equation with both Dirichlet
#and Neumann conditions
#############################################

from dolfin import *
import numpy as np
from tabulate import tabulate
#http://www.math.rutgers.edu/~falk/math575/Boundary-conditions.html


class Poission():
    def __init__(self, h):
        self.y = np.zeros(len(h)); self.y1 = np.zeros(len(h))
        self.x = np.zeros(len(h)); self.h_list = h
        self.L2list = []; self.H1list = []
        self.alpha = 0; self.beta = 0
        self.count = 0

    def set_mesh(self,i):
        self.h = i
        self.mesh = UnitSquareMesh(i, i)

    def calc(self, i, my, output, upwind, imp_norm):
        mesh = self.mesh

        #Defining spaces and functions
        V = FunctionSpace(mesh, 'CG', i)
        u = TrialFunction(V)
        v = TestFunction(V)

        class Left(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], 0)

        class Right(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and near(x[0], 1)

        left = Left(); right = Right()
        #Setting boundary values
        boundaries = FacetFunction("size_t", mesh)
        boundaries.set_all(0)
        left.mark(boundaries,1)
        right.mark(boundaries, 2)
        bc0 = DirichletBC(V, 0, left)
        bc1 = DirichletBC(V, 1, right)
        bcs = [bc0, bc1]

        #Defining and solving variational problem
        u_e = interpolate(Expression('1./(exp(1./my)-my ) * (exp(x[0]/my) - my)', my = my), V)
        f = Constant(0)
        if upwind == True:
            beta_val = 0.5
            beta = Constant(beta_val)
            v = v + beta*v.dx(0)
            a = my * inner(grad(u), grad(v))*dx + u.dx(0)*v*dx #Standard Galerkin
            L = f*v*dx
        else:
            a = my * inner(grad(u), grad(v))*dx + u.dx(0)*v*dx
            L = f*v*dx

        u_ = Function(V)
        solve(a == L, u_, bcs)

        #Norms of the error
        L2 = errornorm(u_e, u_, norm_type='L2', degree_rise = 3)
        H1 = errornorm(u_e, u_, norm_type='H1', degree_rise = 3)

        self.L2list.append(str(L2))
        self.H1list.append(str(H1))

        #plot(u_); interactive()

        if output == True:
            print "----------------------------------"
            print "For %d points and my = %d" % (self.h, my)
            print "L2 Norm = %.5f -----  H1 Norm = %.5f" % (L2, H1)
            print
        if my == 1:

            d = mesh.coordinates()
            self.x[self.count] = np.log(1./self.h)
            self.y[self.count] = np.log( L2) #/ norm(u_e, "H1") )
            self.y1[self.count] = np.log( H1 )
            if imp_norm == True:
                e = u_e-u_
                e = project(e, V)
                i_norm = np.sqrt(mesh.hmin()*norm(e, 'l2')**2 )

            self.count += 1

    def l_square(self, norm, fig):
        A = np.zeros((2, 2))
        b = np.zeros(2)
        if norm == 'H1':
            self.y = self.y1

        A[0][0] = len(self.h_list)
        A[0][1] = np.sum(self.x); A[1][0] = A[0][1]
        A[1][1] = np.sum(self.x*self.x)
        b[0] = np.sum(self.y)
        b[1] = np.sum(self.y*self.x)

        a, b = np.linalg.solve(A, b)
        self.beta = a ; self.alpha = b

        print '####################################################################################\n'
        print '#-------------------------------- Linear Approximation ------------------------------#\n'
        print '                              my = %d' % my[0]
        print '                              alpha = %.4f, beta = %.4f \n' % (prob.alpha, prob.beta)
        print '####################################################################################\n'

        if fig == True:
            import matplotlib.pyplot as plt
            plt.figure(1)
            plt.plot(self.x, b*self.x + a, label='Linear approximation')
            plt.plot(self.x, self.y, 'o', label='Points to be approximated')
            plt.legend(loc = 'upper left')
            plt.show()

    def make_list(self, h):

        k_1 = ['my = 1']; k_10 = ['my = 0.1']; k_100 = ['my = 0.01']
        k_1000 = ['my = 0.001']; k_10000 = ['my = 0.0001']

        for i in range(0, len(self.L2list)-4, 5 ):
            k_1.append(str(self.L2list[i]) )
            k_10.append( str(self.L2list[i+1]) )
            k_100.append( str(self.L2list[i+2]) )
            k_1000.append( str(self.L2list[i+3]) )
            k_10000.append( str(self.L2list[i+4]) )

        table = [k_1, k_10, k_100, k_1000, k_10000]
        headers = ['Values of N']
        for i in h:
            headers.append(str(i))

        print '#------------------------------------ L2 Norm ------------------------------------#\n'
        print tabulate(table, headers, tablefmt="fancy_grid")

        l_1 = ['my = 1']; l_10 = ['my = 0.1']; l_100 = ['my = 0.01']
        l_1000 = ['my = 0.001']; l_10000 = ['my = 0.0001']

        for i in range(0, len(self.H1list)-4, 5 ):
            l_1.append(str(self.H1list[i]) )
            l_10.append( str(self.H1list[i+1]) )
            l_100.append( str(self.H1list[i+2]) )
            l_1000.append( str(self.H1list[i+3]) )
            l_10000.append( str(self.H1list[i+4]) )
        table = [l_1, l_10, l_100, l_1000, l_10000]
        print
        print '#------------------------------------ H1 Norm ------------------------------------#\n'
        print tabulate(table, headers, tablefmt="fancy_grid")
        print

        self.L2list = []
        self.H1list = []


set_log_active(False) #Removing all logging
my = [1*10**-i for i in range(5)]
h = [2**(i+3) for i in range(4)] #5
print my
prob = Poission(h)
for j in [1, 2]:
    print '####################################################################################\n'
    print '#-------------------------------- %d degree elements ------------------------------#\n' % j
    print '####################################################################################\n'
    print
    for i in h:
        for m in my:
            prob.set_mesh(i)
            prob.calc(j, m, output = False, upwind = True, imp_norm = True)
    prob.l_square('L2', fig = False)

    prob.make_list(h)
    prob.count = 0
