############################################
#Author: Andreas Slyngstad
#MEK 4250
#EXERCISE 7.6
#Solving Poission Equation for periodic exact solutions
#############################################

from dolfin import *
import numpy as np
from tabulate import tabulate

def poiseuille(N, v_deg, p_deg):
    mesh = UnitSquareMesh(N, N)
    V = VectorFunctionSpace(mesh, 'CG', v_deg)
    Q = FunctionSpace(mesh, 'CG', p_deg)

    VQ = V*Q
    u, p = TrialFunctions(VQ)
    v, q = TestFunctions(VQ)

    u_e = Expression(("sin(pi*x[1])", "cos(pi*x[0])"))
    p_e = Expression("sin(2*pi*x[0])")
    #plot(u_e, interactive=True)

    f = Expression(("pi*pi*sin(pi*x[1]) - 2*pi*cos(2*pi*x[0])", "pi*pi*cos(pi*x[0])"))

    a = inner(grad(u), grad(v))*dx + p*div(v)*dx + div(u)*q*dx
    L = inner(f, v)*dx

    #F = inner(grad(u), grad(v))*dx +p*div(v)*dx + div(v)*q*dx - inner(f,v)*dx

    #Define boundaries and BCS
    noslip = AutoSubDomain(lambda x: "on_boundary" and not (near(x[0], 0) or near(x[0], 1) ))

    u_bc = DirichletBC(VQ.sub(0),u_e , "on_boundary")
    p_bc = DirichletBC(VQ.sub(1), p_e, "on_boundary")

    bcs = [u_bc, p_bc]

    up = Function(VQ)
    #A = assemble(a); b = assemble(L)
    #solve(A, up.vector(), b, 'gmres')
    solve(a == L, up, bcs)

    u_h, p_h = up.split()
    #plot(u, interactive=True)
    H1 = errornorm(u_e, u_h, norm_type='h1', degree_rise=1)
    L2 = errornorm(p_e, p_h, norm_type='l2', degree_rise=1)
    E.append(H1); h.append(mesh.hmin())
    E_p.append(L2)

    #left = errornorm(u_e, u_h, norm_type='h1', degree_rise=1) + \
    #    errornorm(p_e, p_h, norm_type='l2', degree_rise=1)
    #Left.append(left)
    #right = mesh.hmin()**3*errornorm(u_e, u_h, norm_type='h2', degree_rise=1)
    #te = derivative(u_e, du)

#v_d = [3] ; p_d = [1]
v_d = [4,4,3,3]; p_d = [3,2,2,1]

table = []; headers = ['',]
table2 = [] ; headers2 = ['',]
table3 = []
N = [2**i for i in range(3, 7)]

for i in range(len(v_d)):
    E = []; E_p = []; h = []
    Right = []; Left = []
    print "Solving poiseuille flow for P%d-P%d elements" % (v_d[i], p_d[i])
    for j in N:
        poiseuille(N = j, v_deg = v_d[i], p_deg = p_d[i])
        if i == len(v_d)-1:
            headers.append(str(j))

    li = [];
    li.append(str(v_d[i]) + "-" + str(p_d[i]))
    for k in range(len(E)):
        li.append(E[k])
    table.append(li)

    u_h1 = []; p_l2 = []
    u_h1.append(str(v_d[i]) + "-" + str(p_d[i]))
    p_l2.append(str(v_d[i]) + "-" + str(p_d[i]))
    for k in range(len(E)-1):
        E_u = abs(E[k+1]/E[k]); h_ = abs(h[k+1]/h[k])
        r = np.log(E_u)/np.log(h_)
        u_h1.append(r)
        Ep = abs(E_p[k+1]/E_p[k]);
        r_p = np.log(Ep)/np.log(h_)
        p_l2.append(r_p)
        #print "Convergence from N = %d to N = %d" % (N[k], N[k+1])
        if i == len(v_d)-1:
            headers2.append("Conv %d to %d" % (N[k], N[k+1]))
    table2.append(u_h1)
    table3.append(p_l2)

print tabulate(table, headers)
print "CONVERGENCE RATE FOR VELOCITY"
print tabulate(table2, headers2)
print "CONVERGENCE RATE FOR PRESSURE"
print tabulate(table3, headers2)
