from dolfin import *
import numpy as np
from tabulate import tabulate

def singlespace(N, v_deg, mu, lamb):
    mesh = UnitSquareMesh(N, N)
    V = VectorFunctionSpace(mesh, 'Lagrange', v_deg)

    u = TrialFunction(V)
    v = TestFunction(V)

    mu = mu;
    lamb = lamb;

    f = Expression(("-mu*(-pow(pi,3)*pow(x[0],3)*cos(pi*x[0]*x[1]) \
        -pi*pi*x[1]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))",\
        "-mu*(pow(pi,3)*pow(x[1],3)*cos(pi*x[0]*x[1]) \
            +pi*pi*x[0]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))"), mu = mu)


    analytical = Expression( ("pi*x[0]*cos(pi*x[0]*x[1])", "-pi*x[1]*cos(pi*x[0]*x[1])") )
    u_e = interpolate(analytical, V)

    bc = DirichletBC(V, analytical, "on_boundary")

    #a = mu*inner(grad(u), grad(v))*dx + lamb*inner(div(u), div(v))*dx
    a = mu*inner(grad(u), grad(v))*dx + lamb*inner(div(u), div(v))*dx
    L = dot(f,v)*dx

    u_ = Function(V)
    solve(a == L, u_, bc)


    #plot(u_, interactive=True)
    #plot(u_e, interactive=True)

    L2 = errornorm(u_e, u_, norm_type='l2', degree_rise=3)
    print "The L2 norm of the numerical experiment is %g" % L2
    E.append(L2)
    if lamb == 100 and N == 32:
        ufile_pvd = File("vel_cal.pvd")
        ufile_pvd << u_
        uefile_pvd = File("vel_ex.pvd")
        uefile_pvd << u_e

def mixedspace(N, v_deg, p_deg, mu, lamb):
    mesh = UnitSquareMesh(N, N)
    V = VectorFunctionSpace(mesh, 'Lagrange', v_deg)
    Q = FunctionSpace(mesh, 'Lagrange', p_deg)
    VQ = V*Q

    u, p = TrialFunctions(VQ)
    v, q = TestFunctions(VQ)
    #up = Function(VQ); u, p = split(up)
    #vq = TestFunction(VQ); v, q = split(vq)

    f = Expression(("-mu*(-pow(pi,3)*pow(x[0],3)*cos(pi*x[0]*x[1]) \
        -pi*pi*x[1]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))",\
        "-mu*(pow(pi,3)*pow(x[1],3)*cos(pi*x[0]*x[1]) \
            +pi*pi*x[0]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))"), mu = mu)


    a1 = mu*inner(grad(u), grad(v))*dx + div(v)*p*dx - inner(f,v)*dx
    a2 = lamb*div(u)*q*dx  - p*q*dx
    #L = inner(f,v)*dx
    F = a1 + a2

    analytical = Expression( ("pi*x[0]*cos(pi*x[0]*x[1])", "-pi*x[1]*cos(pi*x[0]*x[1])") )
    bc = DirichletBC(VQ.sub(0), analytical, "on_boundary")

    up_ = Function(VQ)

    solve(lhs(F) == rhs(F), up_, bc)

    u_, p_ = up_.split(True)

    u_e = interpolate(analytical, V)

    #plot(u_, interactive=False)
    #plot(u_e, interactive=False)

    L2 = errornorm(u_e, u_, norm_type='l2', degree_rise=3)
    print "The L2 norm of the numerical experiment is %.5f" % L2
    E.append(L2)


N = [2**i for i in range(3, 7)]
la = [1, 10, 100, 1000, 10000]

h = [1./i for i in N]
table = []; headers = ['N']
tab_con = []; head2 = ['']

solver = "mixed"
for i in range(len(la)):
    E = [];
    for j in N:
        print "N = %d, lamb = %d" % (j, la[i])
        if i == len(la)-1:
            headers.append(str(j))
        if solver == "single":
            singlespace(N = j, v_deg = 2, mu = 1., lamb = la[i])
        if solver == "mixed":
            mixedspace(N = j, v_deg = 2, p_deg = 1, mu = 1., lamb = la[i])



    #L2 norm error
    li = [];
    li.append(str(la[i]) )
    for k in range(len(E)):
        li.append(E[k])
    table.append(li)

    #Convergence
    u_h = []
    u_h.append(str(la[i]))
    for k in range(len(E)-1):
        #Velocity
        E_u = abs(E[k+1]/E[k]); h_ = abs(h[k+1]/h[k])
        r = np.log(E_u)/np.log(h_)
        u_h.append(r)

        if i == len(la)-1:
            head2.append("Conv %d to %d" % (N[k], N[k+1]))

    tab_con.append(u_h)
print tabulate(table, headers)
print tabulate(tab_con, head2)
