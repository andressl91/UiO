from dolfin import *
import numpy as np
import matplotlib.pyplot as plt
#REMARK! Couldn't use higher order space for
#analytical solution due to different array length
#for computing norms

def runcalc(n, k):
    #Spaces and Functions
    mesh = IntervalMesh(n, 0 , 1)
    V = FunctionSpace(mesh, 'CG', 1)
    VH = FunctionSpace(mesh, 'CG', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    f = Expression("k*k*sin(k*pi*x[0])", k = k)

    a = inner(grad(u), grad(v))*dx
    L = f*v*dx


    tol = 1E-8
    left = DirichletBC(V, Constant(0),
            "std::abs(x[0]) < 1e-8")
    right = DirichletBC(V, Constant(0),
            "std::abs(x[0]) > (1-1e-8)")
    bcs = [left, right]
    u_ = Function(V)
    solve(a == L, u_, bcs)



    u_a = interpolate(Expression("sin(k*pi*x[0])", k = k), VH)

    L1 = np.sum(np.abs(u_.vector().array()-u_a.vector().array() ) )
    L2 = norm(u_a, 'l2')

    xeta = Function(V)
    xeta = interpolate(Expression("x[0]"), V)
    xpoint  = xeta.vector().array()
    h_1 = np.abs((xpoint[1]-xpoint[0]))
    C = L1/(h_1*L2)
    print L1
    return C


N = [100, 1000, 10000]
K = [1, 10, 100]
for n in N:
    for k in K:
        print "            N = %d ---- k = %d                    " % (n, k)
        print "                 C = %.3f                          " % runcalc(n, k)
