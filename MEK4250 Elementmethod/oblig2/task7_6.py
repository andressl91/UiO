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

    #WALL SHEAR STRESS
    sides = DomainBoundary() #AutoSubDomain(lambda x: "on_boundary" and near(x[1],0))
    boundaries = FacetFunction("size_t", mesh)
    boundaries.set_all(0)
    sides.mark(boundaries, 1)
    ds = Measure("ds", subdomain_data=boundaries)

    def sigma (u_o):
        return 2.0*mu*sym(grad(u_o)) #- p_o*Identity(len(u_o))

    def eps(u):
        return sym(grad(u))


    #SHEAR STRESS
    n = FacetNormal(mesh)
    t = as_vector((n[1], -n[0]))

    u_e = project(u_e, V)
    #comp = assemble( dot(dot(sigma(u_h), n), t) * ds(1) )
    #exact =  assemble( dot(dot(sigma(u_e), n), t) * ds(1) )

    comp = assemble( dot(sigma(u_h), n)[0] * ds(1) )
    exact =  assemble( dot(sigma(u_e), n)[0] * ds(1) )
    Err = exact - comp
    L2_shear = np.sqrt(Err*Err)
    E_shear.append(L2_shear)


set_log_active(False)

v_d = [4,4,3,3]; p_d = [3,2,2,1]
table = []; headers = ['',]
table2 = [] ; headers2 = ['',]
table3 = [] ; table4 = []
N = [2**i for i in range(2, 7)]
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

        #PLOT ERROR AGAINST h, loglog
        plt.figure(1);
        plt.title("Velocity norm H1")
        plt.loglog(h, E, marker='o', linestyle='--', label = 'Velocity, P%d-P%d ' % ( v_d[i], p_d[i]))
        plt.legend()
        plt.figure(2)
        plt.title("Pressure norm L2")
        plt.loglog(h, E_p, marker='o', linestyle='--', label = 'Pressure, P%d-P%d ' % ( v_d[i], p_d[i]))
        plt.legend()
        plt.figure(3)
        plt.title("Shear STRESS")
        plt.loglog(h, E_shear, marker='o', linestyle='--', label = 'Shear, P%d-P%d ' % ( v_d[i], p_d[i]))
        plt.legend()
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
