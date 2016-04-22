from dolfin import *
import numpy as np

def singlespace(N, count):
    mesh = UnitSquareMesh(N, N)
    V = VectorFunctionSpace(mesh, 'Lagrange', 2)

    u = TrialFunction(V)
    v = TestFunction(V)

    mu = 1.; lamb = 100.;

    f = Expression(("-mu*(-pow(pi,3)*pow(x[0],3)*cos(pi*x[0]*x[1]) \
        -pi*pi*x[1]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))",\
        "-mu*(pow(pi,3)*pow(x[1],3)*cos(pi*x[0]*x[1]) \
            +pi*pi*x[0]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))"), mu = mu)


    analytical = Expression( ("pi*x[0]*cos(pi*x[0]*x[1])", "-pi*x[1]*cos(pi*x[0]*x[1])") )
    u_e = interpolate(analytical, V)

    bc = DirichletBC(V, analytical, "on_boundary")

    a = mu*inner(grad(u), grad(v))*dx + lamb*inner(div(u), div(v))*dx
    L = dot(f,v)*dx

    u_ = Function(V)
    solve(a == L, u_, bc)


    #plot(u_, interactive=True)
    #plot(u_e, interactive=True)

    L2 = errornorm(u_, u_e, norm_type='l2')
    print "The L2 norm of the numerical experiment is %.5f" % L2
    h_single[count] = mesh.hmin()
    E_single[count] = L2

def mixedspace(N, count):
    mesh = UnitSquareMesh(N, N)
    V = VectorFunctionSpace(mesh, 'Lagrange', 2)
    Q = FunctionSpace(mesh, 'Lagrange', 1)
    VQ = V*Q

    up = Function(VQ)
    u, p = split(up)
    vq = TestFunction(VQ)
    v, q = split(vq)

    mu = 1.; lamb = 1.;

    f = Expression(("-mu*(-pow(pi,3)*pow(x[0],3)*cos(pi*x[0]*x[1]) \
        -pi*pi*x[1]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))",\
        "-mu*(pow(pi,3)*pow(x[1],3)*cos(pi*x[0]*x[1]) \
            +pi*pi*x[0]*(2*sin(pi*x[0]*x[1]) + pi*x[0]*x[1]*cos(pi*x[0]*x[1])))"), mu = mu)


    a1 = mu*inner(grad(u), grad(v))*dx + lamb*inner(p, div(v))*dx - inner(f,v)*dx
    a2 = lamb*p*q*dx -lamb*div(u)*q*dx
    F = a1 - a2

    analytical = Expression( ("pi*x[0]*cos(pi*x[0]*x[1])", "-pi*x[1]*cos(pi*x[0]*x[1])") )


    bc = DirichletBC(VQ.sub(0), analytical, "on_boundary")


    up_ = Function(VQ)

    solve(F == 0, up, bc)

    u_, p_ = split(up)
    u_e = project(analytical, V)
    u_ = project(u_, V)

    plot(u_, interactive=False)
    plot(u_e, interactive=False)

    L2 = errornorm(u_, u_e, norm_type='l2')
    print "The L2 norm of the numerical experiment is %.5f" % L2
    h_mixed[count] = mesh.hmin()
    E_mixed[count] = L2


N = [2**i for i in range(3, 5)]

E_single = np.zeros(len(N)); E_mixed = np.zeros(len(N))
h_single = np.zeros(len(N)); h_mixed = np.zeros(len(N))
count = 0
for i in N:
    singlespace(i, count)
    mixedspace(i, count)
    count += 1
print E_single, h_single

r_mixed = np.zeros(len(N)-1); r_single = np.zeros(len(N)-1)

for j in range(len(E_single)-1):
    r_single[j] = np.log(E_single[j+1]-E_single[j]) / np.log(h_single[j+1]-h_single[j])
    r_mixed[j] = np.log(E_mixed[j+1]-E_mixed[j]) / np.log(h_mixed[j+1]-h_mixed[j])

print r_single
