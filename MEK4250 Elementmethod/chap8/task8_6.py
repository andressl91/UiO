from dolfin import *
import time
import numpy as np
import matplotlib.pyplot as plt

def solving_time(A, b, solv, V, tab):

    bc1 = DirichletBC(V, 0, "on_boundary")
    bcs = [bc1]
    [bc.apply(A, b) for bc in bcs]

    u_ = Function(V)
    pc = PETScPreconditioner(solv[1])
    sol = PETScKrylovSolver(solv[0], pc)

    t0 = time.time()
    it = sol.solve(A, u_.vector(), b)

    t1 = time.time()
    tab.append(it)
    return t1-t0


def poission(N, dim, degree):
    if dim == 1:
        mesh = IntervalMesh(N, 0, 1)
    if dim == 2:
        mesh = UnitSquareMesh(N, N)

    V = FunctionSpace(mesh, 'Lagrange', degree)
    u = TrialFunction(V)
    v = TestFunction(V)

    f = Constant(-2)
    #OTHER PROBLEM EXERCISE 8.6
    #u - u'' = f
    #a = u*v*dx + inner(grad(u), grad(v))*dx

    #Poission
    a = inner(grad(u), grad(v))*dx
    L = f*v*dx


    A = assemble(a)
    b = assemble(L)

    u_ = Function(V)

    t0 = time.time()
    it = solve(A, u_.vector(), b, "lu")
    t1 = time.time()
    lu.append(t1-t0)
    lu_table.append(it)

    cg.append(solving_time(A, b, ["cg", "none"], V, cg_table))

    cg_ilu.append(solving_time(A, b, ["cg", "ilu"], V, cgilu_table))

    cg_amg.append(solving_time(A, b, ["cg", "amg"], V, cgamg_table))

    cg_jac.append(solving_time(A, b, ["cg", "jacobi"], V, cgjac_table))

    dof.append(V.dim())



N = [32, 64, 128, 256, 512]

from tabulate import tabulate
set_log_active(False)
for i in [2]:
    lu_table = ["lu"]
    cg_table = ["cg"]; cgilu_table = ["cg_ilu"];
    cgamg_table = ["cg_amg"]; cgjac_table = ["cg_jac"]

    dof = []
    lu = ["lu"];
    cg = ["cg"]; cg_ilu = ["cg_ilu"]; cg_amg = ["cg_amg"]; cg_jac = ["cg_jac"]
    for n in N:
        print "Calculating for n = %d, P-%d elements" % (n, i)
        poission(N = n, dim = i, degree = 1)


    table = [lu, cg,cg_ilu,cg_amg, cg_jac]
    print "TIME vs N, P-%d elements" % i
    print tabulate(table, headers = ['N', 32, 64, 128, 256, 512], tablefmt="fancy_grid")

    print "ITERATIONS vs N, P-%d elements" % i
    table = [lu_table, cg_table,cgilu_table,cgamg_table, cgjac_table]
    print tabulate(table, headers = ['N', 32, 64, 128, 256, 512], tablefmt="fancy_grid")

    np.asarray(lu_table); np.asarray(cg_table)
    np.asarray(cgilu_table); np.asarray(cgamg_table)
    np.asarray(cgjac_table)

    np.asarray(lu); np.asarray(cg)
    np.asarray(cg_ilu); np.asarray(cg_amg)
    np.asarray(cg_jac)

    plt.figure(i)
    plt.plot(dof, lu[1:],label="lu" )
    plt.plot(dof, cg[1:], label="CG")
    plt.plot(dof, cg_ilu[1:], label="CG_ILU")
    plt.plot(dof, cg_amg[1:], label="CG_AMG")
    plt.plot(dof, cg_jac[1:], label="CG_JAC")
    plt.xlabel("DEGREES OF FREEDOM")
    plt.ylabel("TIME ")
    plt.legend(loc=2)
    plt.title("COMPARISON ITERATION TIME,  %d-D problem" % i)

    """
    plt.figure(i+10)
    plt.plot(N, lu_table[1:], label="lu" )
    plt.plot(N, cg_table[1:], label="cg" )
    plt.plot(N, cgilu_table[1:], label="cg_ilu" )
    plt.plot(N, cgamg_table[1:], label="cg_amg" )
    plt.plot(N, cgjac_table[1:], label="cg_jac" )
    plt.xlabel("DEGREES OF FREEDOM")
    plt.ylabel("TIME ")
    plt.legend(loc=2)
    plt.title("ITERATIONS  %d-D problem" % i)
    """

plt.show()
