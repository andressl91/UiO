from dolfin import *
from numpy import log
import matplotlib.pyplot as plt

def condiff(N, mu, i, last, typ):
    mesh = IntervalMesh(N, 0, 1)

    V = FunctionSpace(mesh, 'Lagrange', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    f = Constant(0)

    if typ == "artificial":
        h_dist = 1./N
        beta = Constant(0.5)
        a = -u.dx(0)*v*dx + (mu +beta*h_dist)*u.dx(0)*v.dx(0)*dx
        L = f*v*dx

    if typ == "naive":
        a = -u.dx(0)*v*dx + mu*inner(grad(u),grad(v))*dx
        L = f*v*dx

    if typ == "petgal":
        beta = Constant(0.5)
        v = v + beta*v.dx(0)
        a = -u.dx(0)*v*dx + mu*inner(grad(u),grad(v))*dx
        L = f*v*dx



    u_e = Expression("(exp(-x[0]/mu) - 1 ) / (exp(-1/mu) - 1 ) ", mu = mu)
    u_e = project(u_e, V)

    bc0 = DirichletBC(V, 0, "near(x[0], 0)")
    bc1 = DirichletBC(V, 1, "near(x[0], 1)")
    bcs = [bc0, bc1]

    u_ = Function(V)

    solve(a == L, u_, bcs)
    h.append(mesh.hmin())
    E.append(errornorm(u_e, u_, norm_type='l2', degree_rise = 3)  )

    u_ = u_.vector().array()
    u_e = u_e.vector().array()
    x_val = mesh.coordinates()


    plt.figure(1)
    #plt.subplot("%d" %(220 + i))
    plt.title("Con-Vec")
    plt.plot(u_[::-1], label='mu=%g' % mu)
    #plt.plot(x_val, u_e[::-1], label='Exact')
    if i == last:
        u_o = Expression("1")
        u_o = project(u_o, V)
        u_o = u_o.vector().array()
        plt.plot(u_o[::-1],  label = 'mu=%g' % 0)
        #plt.axis([0, 1, 0, 1.2])
        plt.legend(loc=4)
        plt.show()


N = 100
#mu = [1, 0.1, 0.01, 0.001]
mu = [1, 0.5, 0.1, 0.05, 0.01]
#N = [2**i for i in range(3, 13)]
#mu = [0.01]

set_log_active(False)
for k in ["artificial"]:
    E = []; h = []
    for i in range(len(mu)):
        condiff(N, mu[i], i, len(mu)-1, typ = k)
    for j in range(len(E)-1):
        r = log(E[j+1]/E[j]) / log(h[j+1]/h[j])
        print r
