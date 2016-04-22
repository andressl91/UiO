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
    def __init__(self, h, lamb):
        self.y = np.zeros(len(h)); self.y1 = np.zeros(len(h))
        self.x = np.zeros(len(h)); self.h_list = h
        self.L2list = []; self.H1list = []
        self.alpha = 0; self.beta = 0
        self.count = 0

    def set_mesh(self,i):
        self.h = i
        self.mesh = UnitSquareMesh(i, i)

    def calc(self, i, lamb, space, output=True):
        mesh = self.mesh
        mu = 1.;

        f = Expression(("-mu*(-pow(pi,3)*pow(x[0],3)*cos(pi*x[0]*x[1]) \
            -pi*pi*x[1]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))",\
            "-mu*(pow(pi,3)*pow(x[1],3)*cos(pi*x[0]*x[1]) \
                +pi*pi*x[0]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))"), mu = mu)

        if space == "singelspace":
            #Defining spaces and functions
            V = VectorFunctionSpace(mesh, 'Lagrange', i)
            V1 = VectorFunctionSpace(mesh, 'Lagrange', i+1)
            u = TrialFunction(V)
            v = TestFunction(V)

            analytical = Expression( ("pi*x[0]*cos(pi*x[0]*x[1])", "-pi*x[1]*cos(pi*x[0]*x[1])") )
            u_e = interpolate(analytical, V1)

            bc = DirichletBC(V, analytical, "on_boundary")

            a = mu*inner(grad(u), grad(v))*dx + lamb*inner(div(u), div(v))*dx
            L = dot(f,v)*dx

            u_ = Function(V)
            solve(a == L, u_, bc)

            L2 = errornorm(u_e, u_, norm_type='L2', degree_rise = 3)
            H1 = errornorm(u_e, u_, norm_type='H1', degree_rise = 3)

            if output == True:
                print "----------------------------------"
                print "Velocity P%s" % str(i)
                print "For %d points and lambda = %d" % (self.h, lamb)
                print "L2 Norm = %.5f -----  H1 Norm = %.5f" % (L2, H1)
                print

        else:
            V = VectorFunctionSpace(mesh, 'Lagrange', i+1)
            V1 = VectorFunctionSpace(mesh, 'Lagrange', i+2)
            Q = FunctionSpace(mesh, 'Lagrange', i)
            VQ = V*Q

            up = Function(VQ)
            u, p = split(up)
            vq = TestFunction(VQ)
            v, q = split(vq)

            a1 = mu*inner(grad(u), grad(v))*dx + lamb*inner(p, div(v))*dx - inner(f,v)*dx
            a2 = lamb*p*q*dx -lamb*div(u)*q*dx
            F = a1 - a2

            analytical = Expression( ("pi*x[0]*cos(pi*x[0]*x[1])", "-pi*x[1]*cos(pi*x[0]*x[1])") )
            u_e = interpolate(analytical, V1)

            bc = DirichletBC(VQ.sub(0), analytical, "on_boundary")


            up_ = Function(VQ)

            solve(F == 0, up, bc)

            u_, p_ = split(up)
            u_ = project(u_, V)

            L2 = errornorm(u_e, u_, norm_type='L2', degree_rise = 3)
            H1 = errornorm(u_e, u_, norm_type='H1', degree_rise = 3)

            if output == True:
                print "----------------------------------"
                print "Velocity P%s, Pressure P%s " % (str(i+1), str(i))
                print "For %d points and lambda = %d" % (self.h, lamb)
                print "L2 Norm = %.5f -----  H1 Norm = %.5f" % (L2, H1)
                print

        #Norms of the error

        self.L2list.append(str(round(L2,5) ))
        self.H1list.append(str(round(H1,5) ))

        self.x[self.count] = np.log(1./self.h)
        self.y[self.count] = np.log( L2 )
        self.y1[self.count] = np.log( H1 )
        self.count += 1


    def l_square(self, lamb, norm,fig):
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
        print '                       Norm = %s     lambda = %d' % (norm, lamb)
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

        k_1 = ['lambda = 1']; k_10 = ['lambda = 10']; k_100 = ['lambda = 100']
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

        l_1 = ['lambda = 1']; l_10 = ['lambda = 10']; l_100 = ['lambda = 1']

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

h = [2**(i+3) for i in range(3)]
lamb = [1, 10, 100]
prob = Poission(h, lamb)
for s in ["singlespace", "mixedspace"]:
    print '###########################################################\n'
    print '#-------------------- %s_space -----------------#\n' % s
    print '###########################################################\n'
    for j in [1, 2]:
        print '###########################################################\n'
        print '#-------------------- %d degree elements -----------------#\n' % j
        print '###########################################################\n'
        print
        for l in lamb:
            for i in h:
                prob.set_mesh(i)
                prob.calc(j, l, s, output = True)

        print '###########################################################\n'
        print '#-------------------- Linear Approximation ---------------#\n'
        for g in ['L2', 'H1']:
            prob.l_square(l, g, fig = False)
        print
        print '###########################################################\n'
        #prob.make_list(h)
        prob.count = 0


"""
h = [2**(i+3) for i in range(3)]
lamb = [1, 10, 100]
prob = Poission(h, lamb)
for s in ["singlespace", "mixedspace"]:
    print '###########################################################\n'
    print '#-------------------- %s_space -----------------#\n' % s
    print '###########################################################\n'
    for j in [1, 2]:
        print '###########################################################\n'
        print '#-------------------- %d degree elements -----------------#\n' % j
        print '###########################################################\n'
        print
        for i in h:
            for l in lamb:
                prob.set_mesh(i)
                prob.calc(j, l, s, output = False)

        print '###########################################################\n'
        print '#-------------------- Linear Approximation ---------------#\n'
        for l in ['L2', 'H1']:
            prob.l_square(l, fig = False)
        print
        print '###########################################################\n'
        #prob.make_list(h)
        prob.count = 0
"""
