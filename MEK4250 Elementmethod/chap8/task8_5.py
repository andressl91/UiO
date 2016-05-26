from dolfin import *
import numpy as np
import matplotlib.pyplot as plt

def exact(N, A, B):
    mesh = IntervalMesh(N-1, 0, 1)
    V = FunctionSpace(mesh, 'Lagrange', 1)
    u = TrialFunction(V)
    v = TestFunction(V)

    f = Constant(2)

    a = inner(grad(u), grad(v))*dx
    L = f*v*dx

    u_ = Function(V)
    bc1 = DirichletBC(V, A, "std::abs(x[0]) < DOLFIN_EPS")
    bc2 = DirichletBC(V, B, "std::abs(x[0]) > 1 - DOLFIN_EPS")
    bcs = [bc1, bc2]

    A,b = assemble_system(a,L, bcs)
    solve(A,u_.vector(), b)
    """
    print "FENICS A-MATRIX"
    print A.array()
    print "FENICS B-MATRIX"
    print b.array()
    """

    return u_.vector().array()[::-1]

def poission(N, S, T):

    h = 1./(N-1)
    A = np.zeros([N,N])

    for i in range(1, N-1):
        A[i][i] = 2/h
        A[i][i-1] = -1/h
        A[i][i+1] = -1/h
        if i == 1:
            A[i-1][i-1] = 1
            A[i][i+1] = -1/h
            A[i][i-1] = 0
        if i == N-2:
            A[i+1][i+1] = 1
            A[i][i+1] = 0

    b = np.zeros(N)
    for i in range(0, N-1):
        b[i] = 2*h
        if i == 1:
            b[i-1] = S;
            b[i] = 2*h + S/h
        if i == N-1:
            b[i-1] = 2*h + T/h
            b[i] = T



    """
    print "OWN A-MATRIX"
    print A

    print "OWN b-MATRIX"
    print b
    """
    D = np.diag(np.diag(A))
    R = A - D
    D_1 = np.linalg.inv(D)
    u_prev = np.random.randint(5, size=N)

    eigenvalues = np.sort(np.linalg.eigvals(A))

    # compute eigenvalues and tau
    lambda_max, lambda_min = eigenvalues[-1], eigenvalues[0]
    print "lambda_max ", lambda_max, " lambda_min ", lambda_min
    tau = 2/(lambda_max + lambda_min)
    kappa = lambda_max/lambda_min

    ######
    estimate = np.log(10**6)/(np.log( (kappa-1)/(kappa+1) ) )
    print "ITERATION ESTIMATE %f" % estimate


    num_it = 1

    u = D_1.dot(b) - D_1.dot(R.dot(u_prev))
    r = A.dot(u) - b
    residue = np.sqrt((r.dot(r.T)))
    #limit = 10**4*residue
    #print limit, residue
    residue = 1

    while 1e-6 < residue and num_it < 1000:
        u = D_1.dot(b) - D_1.dot(R.dot(u_prev))
        r = A.dot(u) - b
        residue = np.sqrt((r.dot(r.T)))
        print residue
        u_prev = u
        num_it += 1

    print "NUMBER ITERATION %d" % num_it
    x = np.linspace(0, 1, N)
    u_e = exact(N, S, T)


    plt.figure(1)
    plt.plot(x, u_e, label = "Exact")
    plt.plot(x, u, "*",label = "%d it" % num_it)
    #plt.axis([0, 1, 0, 0.4])
    plt.legend()
    plt.show()

poission(N = 6, S = 0, T = 0)
