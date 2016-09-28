
import numpy as np
from scipy.linalg import hilbert

def make(n):
    def hilb(i, j):
        return 1 / ((i+1) + (j+1) - 1)
    #n = 10
    Pnum = np.eye(n+1)

    Hilbert = np.fromfunction(lambda i, j: 1 / ((i+1) + (j+1) - 1), (n+1, n+1), dtype=float)

    altsign = [(-1)**i for i in range(n+1)]
    altsign = np.asarray(altsign)
    altsignmatrix = np.matrix([altsign*(-1)**i for i in range(n+1)])
    Nnum = Hilbert + np.multiply(altsignmatrix, Hilbert)
    #print Nnum

    for i in range(1, n+1):
        A = Nnum[:i, :i]
        b = (-1*(Nnum[:i, i]))
        x = np.linalg.solve(A, b)
        x = np.transpose(x[:i])
        Pnum[:i, i] = x

    #print "Numerical matrix N"
    #print Nnum
    #print "Numerical matrix P"
    #print Pnum

    #Symbolic
    from sympy import *
    from mpmath import norm
    from IPython.display import display

    def f(i, j):
        if i % 2 == 0:
            return 1 if j % 2 == 0 else -1
        if i % 2 == 1:
            return -1 if j % 2 == 0 else 1

    altsignmatrix = Matrix(n+1, n+1, f)
    mat = Matrix(n+1, n+1, hilb)
    Nnumsym = eye(n+1)
    Pnumsym = eye(n+1)
    for i in range(n+1):
        for j in range(n+1):
            Nnumsym[i, j] = mat[i, j] + altsignmatrix[i,j]*mat[i, j]


    for i in range(1, n+1):
        A = Nnumsym[:i, :i]
        b = (-1*(Nnumsym[:i, i]))
        x = A.LUsolve(b)
        Pnumsym[:i, i] = x
    """
    print "Symbolic matrix N"
    display(Nnumsym)
    print "Symbolic matrix P"
    display(Pnumsym)
    """

    #basis = [lambda x: x**i for i in range(n)]
    #coeff = [i for i in range(n)]
    #coeff = Pnum[:, 3]
    """
    evalPoly = lambda coeff, x: sum((x**power) * coeff for power, coeff in enumerate(coeff))
    import matplotlib.pyplot as plt
    x_points = np.linspace(-1, 1, 500)
    error = []
    plt.figure(1) #Numerical
    plt.title("First %d polynomials \n Numerical approximation" % n)
    for i in range(n):
        plt.plot(evalPoly(Pnum[:,i], x_points))
    plt.figure(2)
    plt.title("First %d polynomials \n Symbolic approximation" % n)
    for i in range(n):
        plt.plot(evalPoly(Pnumsym[:,i], x_points))
    """
    print "FOR N = %f CONDITION NUMBER NUMERIC  %.4f" % (n, np.linalg.cond(Pnum))
    #print Pnumsym.condition_number()
    #print "FOR N = %f CONDITION SYMBOLIC NUMERIC  %g" % (n, Pnumsym.condition_number())
    #print "CONDITION NUMBER SYMBOLIC %.4f" % (np.linalg.cond(Pnumsym))
    """
    plt.figure(3) #Compare symbolic and numerical for given n
    plt.title("Legendre polynomial %d" % (n))
    plt.plot(evalPoly(Pnum[:,n], x_points), label="Numerical")
    plt.plot(evalPoly(Pnumsym[:,n], x_points), label="Symbolic")
    plt.legend()
    plt.show()
    """
for i in range(25):
    make(i)
