from dolfin import *
import matplotlib.pyplot as plt

def runcalc(n, k):
    #Spaces and Functions
    mesh = IntervalMesh(n, 0 , 1)
    V = FunctionSpace(mesh, 'CG', 1)
    VH = FunctionSpace(mesh, 'CG', 2)
    u = TrialFunction(V)
    v = TestFunction(V)

    k = 1
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

    xeta = Function(V)
    xeta = interpolate(Expression("x[0]"), V)
    xpoint  = xeta.vector().array()

    u_ = u_.vector().array()
    return u_, xpoint

N = [100, 1000, 10000]
for i in N:
    u, xpoint = runcalc(i, 1)
    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    ax1.plot(xpoint, u)

ax1.set_xlabel("x")
ax1.set_ylabel("u(x)")
ax1.legend(["N =  %d" % p for p in N])

plt.show()
    #ana_u = interpolate(Expression("sin(k*pi*x[0])", k = k), VH)
    #plot(func); interactive()
