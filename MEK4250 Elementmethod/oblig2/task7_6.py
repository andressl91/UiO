############################################
#Author: Andreas Slyngstad
#MEK 4250
#EXERCISE 7.6
#Solving Poission Equation for periodic exact solutions
#############################################

from dolfin import *
import matplotlib.pyplot as plt
import numpy as np
from tabulate import tabulate

def poiseuille(N, v_deg, p_deg, mu):
    mesh = UnitSquareMesh(N, N)
    V = VectorFunctionSpace(mesh, 'Lagrange', v_deg)
    Q = FunctionSpace(mesh, 'Lagrange', p_deg)

    VQ = V*Q
    u, p = TrialFunctions(VQ)
    v, q = TestFunctions(VQ)

    u_e = Expression(("sin(pi*x[1])", "cos(pi*x[0])"))
    p_e = Expression("sin(2*pi*x[0])")
    u_e = interpolate(u_e,V)

    f = Expression(("pi*pi*sin(pi*x[1]) - 2*pi*cos(2*pi*x[0])", "pi*pi*cos(pi*x[0])"))

    a = inner(grad(u), grad(v))*dx + div(v)*p*dx + div(u)*q*dx
    L = inner(f, v)*dx

    #F = inner(grad(u), grad(v))*dx +p*div(v)*dx + div(v)*q*dx - inner(f,v)*dx

    #Define boundaries and BCS

    u_bc = DirichletBC(VQ.sub(0), u_e , "on_boundary")
    p_bc = DirichletBC(VQ.sub(1), p_e, "on_boundary")

    bcs = [u_bc, p_bc]

    up = Function(VQ)
    #A = assemble(a); b = assemble(L)
    #A, b = assemble_system(a, L, bcs)
    #solve(A, up.vector(), b)
    solve(a == L, up, bcs)

    u_h, p_h = up.split(True)
    #plot(u_h.dx(1))
    #plot(u_e.dx(1))
    #plot(u_e)
    #plot(u_h);
    #interactive()
    H1 = errornorm(u_e, u_h, norm_type='h1', degree_rise=2)
    L2 = errornorm(p_e, p_h, norm_type='l2', degree_rise=2)
    E.append(H1); h.append(mesh.hmin())
    E_p.append(L2)

    #WALL SHEAR STRESS
    sides = AutoSubDomain(lambda x:  ( x[1] < DOLFIN_EPS or x[1] > 1 - DOLFIN_EPS ) )
    class Wall(SubDomain):
        def inside(self, x, on_boundary):
            return  x[1] < DOLFIN_EPS or x[1] > 1 - DOLFIN_EPS

    wall = Wall()
    boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    wall.mark(boundaries, 1)
    #ds = Measure("ds", subdomain_data=boundaries)
    ds = Measure("ds", domain = mesh, subdomain_data=boundaries)

    #V1 = VectorFunctionSpace(mesh, 'Lagrange', v_deg+1)
    u_e = interpolate(u_e, V)
    p_e = interpolate(p_e, Q)

    def sigma (u_o, p_o):
        return 2.0*mu*sym(grad(u_o)) - p_o*Identity(len(u_o))

    def eps(u):
        return sym(grad(u))

    n = FacetNormal(mesh)
    t = as_vector((n[1], -n[0]))

    #comp = assemble( dot(dot(sigma(u_h), n), t) * ds(1) )
    #exact =  assemble( dot(dot(sigma(u_e), n), t) * ds(1) )

    #comp = assemble( dot(t, dot(sigma(u_h, p_h), n)) * ds(1) )
    #exact = assemble( dot(t, dot(sigma(u_e, p_e), n)) * ds(1) )

    #comp = assemble( dot(sigma(u_h, p_h), n)[0] * ds(1) )
    #exact =  assemble( dot(sigma(u_e, p_e), n)[0] * ds(1) )
    #Err = exact - comp
    Err = sqrt( assemble( ( u_h.dx(1)[0]-u_e.dx(1)[0] )**2 * ds(1)) )
    E_shear.append(Err)


set_log_active(False)

v_d = [4,4,3,3]; p_d = [3,2,2,1]
table = []; headers = ['',]
table2 = [] ; headers2 = ['',]
table3 = [] ; table4 = []
N = [2**i for i in range(3, 7)]
for i in range(len(v_d)):
    E = []; E_p = []; h = []; E_shear = []
    Right = []; Left = []
    if MPI.rank(mpi_comm_world()) == 0:
        print "Solving poiseuille flow, P%d-P%d elements" % (v_d[i], p_d[i])
    for j in N:
        if MPI.rank(mpi_comm_world()) == 0:
            print "N = %d" % j
        poiseuille(N = j, v_deg = v_d[i], p_deg = p_d[i], mu = 1.)
        if i == len(v_d)-1:
            headers.append(str(j))

    if MPI.rank(mpi_comm_world()) == 0:
        #L2 norm error
        li = [];
        li.append(str(v_d[i]) + "-" + str(p_d[i]))
        for k in range(len(E)):
            li.append(E[k])
        table.append(li)

        #Convergence u, p, shear
        u_h1 = []; p_l2 = []; u_shearl2 = []
        u_h1.append(str(v_d[i]) + "-" + str(p_d[i]))
        p_l2.append(str(v_d[i]) + "-" + str(p_d[i]))
        u_shearl2.append(str(v_d[i]) + "-" + str(p_d[i]))
        for k in range(len(E)-1):
            #Velocity
            E_u = abs(E[k+1]/E[k]); h_ = abs(h[k+1]/h[k])
            r = np.log(E_u)/np.log(h_)
            u_h1.append(r)
            #Pressure
            Ep = abs(E_p[k+1]/E_p[k]);
            r_p = np.log(Ep)/np.log(h_)
            p_l2.append(r_p)
            #Shear
            E_sh = abs(E_shear[k+1]/E_shear[k]); h_ = abs(h[k+1]/h[k])
            r_sh = np.log(E_sh)/np.log(h_)
            u_shearl2.append(r_sh)

            #print "Convergence from N = %d to N = %d" % (N[k], N[k+1])
            if i == len(v_d)-1:
                headers2.append("Conv %d to %d" % (N[k], N[k+1]))


    E = np.asarray(E); E_p = np.asarray(E_p)
    #PLOT ERROR AGAINST h, loglog
    plt.figure(4);
    #plt.subplot(211)
    plt.title("Velocity norm H1 + pressure norm L2")
    plt.loglog(h, E+E_p, marker='o', linestyle='--', label = 'Elements, P%d-P%d ' % ( v_d[i], p_d[i]))
    plt.legend(loc=4)
    plt.savefig("comb1.png")

    plt.figure(1);
    #plt.subplot(211)
    plt.title("Velocity norm H1")
    plt.loglog(h, E, marker='o', linestyle='--', label = 'Velocity, P%d ' % ( v_d[i]))
    plt.legend(loc=4)
    plt.savefig("velocity1.png")

    plt.figure(2)
    plt.title("Pressure norm L2")
    plt.loglog(h, E_p, marker='o', linestyle='--', label = 'Pressure, P%d' % ( p_d[i]))
    plt.legend(loc=4)
    plt.savefig("pressure1.png")

    plt.figure(3)
    plt.title("Shear STRESS")
    plt.loglog(h, E_shear, marker='o', linestyle='--', label = 'Shear, P%d-P%d ' % ( v_d[i], p_d[i]))
    plt.legend(loc=4)
    plt.savefig("shear.png")

    table2.append(u_h1)
    table3.append(p_l2)
    table4.append(u_shearl2)

plt.show()
if MPI.rank(mpi_comm_world()) == 0:
    print "L2 norm VELOCITY"
    print tabulate(table, headers)
    print "CONVERGENCE RATE FOR VELOCITY"
    print tabulate(table2, headers2)
    print "CONVERGENCE RATE FOR PRESSURE"
    print tabulate(table3, headers2)
    print "CONVERGENCE RATE FOR SHEAR STRESS"
    print tabulate(table4, headers2)
