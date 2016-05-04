############################################
#Author: Andreas Slyngstad
#MEK 4250
#EXERCISE 1
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

    def calc(self, i, k, l, output=True):
        mesh = self.mesh

        #Defining spaces and functions
        V = FunctionSpace(mesh, 'CG', i)
        u = TrialFunction(V)
        v = TestFunction(V)

        class Dirichlet(SubDomain):
            def inside(self, x, on_boundary):
                return on_boundary and( near(x[0], 0) or near(x[0], 1) )

        diri = Dirichlet()
        #Setting boundary values
        boundaries = FacetFunction("size_t", mesh)
        boundaries.set_all(0)
        diri.mark(boundaries,1)
        bc0 = DirichletBC(V, 0, diri)

        #Defining and solving variational problem
        V_1 = FunctionSpace(mesh, 'CG', i+2)
        u_e = interpolate(Expression('sin(k*pi*x[0])*cos(l*pi*x[1])', k=k, l=l), V_1)
        f = Expression("((pi*pi*k*k)+(pi*pi*l*l))*sin(pi*k*x[0])*cos(pi*l*x[1])",k=k,l=l)
        a = inner(grad(u), grad(v))*dx
        L = f*v*dx

        u_ = Function(V)
        solve(a == L, u_, bc0)

        #Norms of the error
        L2 = errornorm(u_e, u_, norm_type='L2', degree_rise = 3)
        H1 = errornorm(u_e, u_, norm_type='H1', degree_rise = 3)
        self.L2list.append(str(round(L2,4) ))
        self.H1list.append(str(round(H1,4) ))

        if output == True:
            print "----------------------------------"
            print "For %d points and k, l = %d" % (self.h, k)
            print "L2 Norm = %.5f -----  H1 Norm = %.5f" % (L2, H1)
            print
        if k == 1:
            d = mesh.coordinates()
            self.x[self.count] = np.log(1./self.h)
            self.y[self.count] = np.log( L2 )
            self.y1[self.count] = np.log( H1 )
            self.count += 1

    def l_square(self, norm ,fig):
        A = np.zeros((2, 2))
        b = np.zeros(2)

        mid = self.y #hold y values if norm = H1
        test = self.y #Holds L2 errornorms
        if norm == 'H1':
            self.y = self.y1
            test = self.y1 #Holds H! errornorms

        A[0][0] = len(self.h_list)
        A[0][1] = np.sum(self.x); A[1][0] = A[0][1]
        A[1][1] = np.sum(self.x*self.x)
        b[0] = np.sum(self.y)
        b[1] = np.sum(self.y*self.x)

        a, b = np.linalg.solve(A, b)
        self.beta = a ; self.alpha = b

        print
        print '                       Norm = %s     k_l = %d' % (norm ,1)
        print '                     alpha = %.4f, Constant = %.4f \n' % (prob.alpha, exp(prob.beta))
        for i in range(len(test)):
            print 'Errornorm (u-u_h) < C*h^(alpha) is %s for N = %d' %(test[i]<b*self.h**a, self.h_list[i])

        if fig == True:
            import matplotlib.pyplot as plt
            plt.figure(1)
            plt.plot(self.x, b*self.x + a, label='Linear approximation')
            plt.plot(self.x, self.y, 'o', label='Points to be approximated')
            plt.legend(loc = 'upper left')
            plt.show()
        self.y = mid

    def make_list(self, h):

        k_1 = ['k_l = 1']; k_10 = ['k_1 = 10']; k_100 = ['k_l = 100']
        for i in range(0, len(self.L2list)-2, 3 ):
            k_1.append(str(self.L2list[i]) )
            k_10.append( str(self.L2list[i+1]) )
            k_100.append( str(self.L2list[i+2]) )

        table = [k_1, k_10, k_100]
        headers = ['Values of N']
        for i in h:
            headers.append(str(i))

        print '#----------------------- L2 Norm ---------------------#\n'
        print tabulate(table, headers, tablefmt='rst')

        l_1 = ['k_l = 1']; l_10 = ['k_1 = 10']; l_100 = ['k_l = 100']

        for i in range(0, len(self.H1list)-2, 3 ):
            l_1.append(str(self.H1list[i]) )
            l_10.append( str(self.H1list[i+1]) )
            l_100.append( str(self.H1list[i+2]) )
        table = [l_1, l_10, l_100]
        print
        print '#----------------------- H1 Norm ---------------------#\n'
        print tabulate(table, headers, tablefmt='rst')
        print

        self.L2list = []
        self.H1list = []


set_log_active(False) #Removing all logging
kl = [1, 10, 100]
h = [2**(i+3) for i in range(4)]

prob = Poission(h)
for j in [1, 2]:
    print '###########################################################\n'
    print '#-------------------- %d degree elements -----------------#\n' % j
    print '###########################################################\n'
    print
    for i in h:
        for k in kl:
            prob.set_mesh(i)
            prob.calc(j, k, k, output = False)

    print '###########################################################\n'
    print '#-------------------- Linear Approximation ---------------#\n'
    for l in ['L2', 'H1']:
        prob.l_square(l, fig = False)
    print
    print '###########################################################\n'
    prob.make_list(h)
    prob.count = 0
