from dolfin import *


def condiff(N, mu):
    mesh = IntervalMesh(N, 0, 1)

    V = FunctionSpace(mesh, 'Lagrange', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    a = -grad(u)*v*dx + mu*grad(u)*grad(v)*dx
    L = 0
    u_e = Expression("(exp(-x[0]/mu) - 1 ) / (exp(-1/mu) - 1 ) ", mu = mu)
    u_e = project(u_e, V)

    bc0 = DirichletBC(V, 0, near(x[0], 0))
    bc1 = DirichletBC(V, 1, near(x[0], 1))
    bcs = [bc0, bc1]

    u_ = Function(V)

    solve(a == L, u_, bcs)

    plot(u_, interactive=True)


N = 10
mu = 1
